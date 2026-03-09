#include "pm2_bailinbailout.h"

/*<===========================================================================>

    Winter 2026
    Author: Kyryll Kotyk
    Advisor: Munehiro Fukuda

    This file calls all phase implementations using given parameters

    Style Note:
        One parameter per line for function signatures
        Parantheses are closed on a new line, followed by body brace opening
        This rule does not apply if there are no parameters

<===========================================================================>*/

void BailInBailOut::runBenchmarkInAndOut(
    /* Simulation Specifiers */
    const unsigned int runs,
    const unsigned int timesteps,
    const unsigned int logEvery,
    const uint64_t baseSeed,

    /* Structural Parameters*/
    const unsigned int bankCountTotal,
    const unsigned int firmCountTotal,
    const unsigned int workerCountTotal,
    // How many workers each bank has
    const unsigned int bankWorkerCount,

    /* System Parameters */
    const unsigned int initialBankLiquidity,
    const unsigned int initialFirmLiquidity,
    const unsigned int initialProductionCost,
    const unsigned int wage,
    const unsigned int bankEmployeeWage,

    /* Multipliers (in %, 95 = 0.95 multiplier) */
    unsigned short wageConsumptionPercent,
    unsigned short profitMultiplierMin,
    unsigned short profitMultiplierMax,
    unsigned short shockMultiplierMin,
    unsigned short shockMultiplierMax,
    unsigned short minInterestRate,
    unsigned short maxInterestRate,

    /* Banking Interaction Parameters */
    unsigned short interbankDensity,
    const unsigned short maxInterbankLenderSamplingK,
    // Maximum percentage of the bank's liquidity it can offer as a loan
    unsigned short maxInterbankLoanPercent,
    unsigned short maxFirmLoanPercent,
    unsigned short firmLenderDegree,
    unsigned short firmRepayPercent,
    unsigned short bankRepayPercent,

    unsigned int interventionDelay,

    /* Bail-In Parameters */
    unsigned short bailInCoveragePercent,

    /* Bail-Out Parameters */
    unsigned short bailOutCoveragePercent,

    bool debug
) {

    // <!===================================================================!>
    // <!=========================!Error Detection!=========================!>
    // <!===================================================================!>

    if (!errorDetectionAndClamping(runs, timesteps, logEvery, baseSeed,
        bankCountTotal, firmCountTotal, workerCountTotal, initialBankLiquidity,
        initialFirmLiquidity, initialProductionCost, wage,
        wageConsumptionPercent, profitMultiplierMin, profitMultiplierMax,
        shockMultiplierMin, shockMultiplierMax, minInterestRate, maxInterestRate,
        interbankDensity, maxInterbankLenderSamplingK, maxInterbankLoanPercent,
        maxFirmLoanPercent, firmLenderDegree, firmRepayPercent,
        bankRepayPercent, interventionDelay, bailInCoveragePercent, bailOutCoveragePercent)) {
        return;
    }

    // <!===================================================================!>
    // <!=========================!MPI Information!=========================!>
    // <!===================================================================!>

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

    /*if (mpiRank == 0) {
        printf("sizeof(e3InterbankLoanAcceptance) = %zu\n", sizeof(e3InterbankLoanAcceptance));
        printf("sizeof(e3InterbankLoanRequest)    = %zu\n", sizeof(e3InterbankLoanRequest));
        printf("sizeof(e1InterbankRepaymentMessage) = %zu\n", sizeof(e1InterbankRepaymentMessage));
        printf("sizeof(e3NeighborLiquidityMessage) = %zu\n", sizeof(e3NeighborLiquidityMessage));
        printf("sizeof(FirmLoanRequest) = %zu\n", sizeof(FirmLoanRequest));
        printf("sizeof(d2FirmLoanAcceptance) = %zu\n", sizeof(d2FirmLoanAcceptance));
    }*/

    int mpiSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    // <!===================================================================!>
    // <!=========================!Range and Count!=========================!>
    // <!===================================================================!>

    // Banks
    unsigned int bankCountForRank, bankGlobalStartIndex;
    pair<unsigned int, unsigned int> bankPair = computeRange(bankCountTotal,
        mpiRank, mpiSize);
    bankCountForRank = bankPair.first;
    bankGlobalStartIndex = bankPair.second;

    // Firms
    unsigned int firmCountForRank, firmGlobalStartIndex;
    pair<unsigned int, unsigned int> firmPair = computeRange(firmCountTotal,
        mpiRank, mpiSize);
    firmCountForRank = firmPair.first;
    firmGlobalStartIndex = firmPair.second;

    // Workers
    unsigned int workerCountForRank, workerGlobalStartIndex;
    pair<unsigned int, unsigned int> workerPair = computeRange(workerCountTotal,
        mpiRank, mpiSize);
    workerCountForRank = workerPair.first;
    workerGlobalStartIndex = workerPair.second;

    uint64_t totalBailOutMoney = 0, totalBailInMoney = 0, totalGraceMoney = 0;
    // Pre-F snapshots (captured after E.3, before f1f2 runs)
    vector<int64_t> preFBankLiquidity;
    vector<vector<DebtEntry>> preFBankDebts;
    vector<unsigned int> preFBankDistressCount;
    // Phase-level snapshots for mid-step correctness verification
    // snapPostA: after A.5.5 employee wages deducted (firm repayments NOT yet applied)
    // snapPostC: after c1SendMoneyToBanks (deposits + repayments received)
    // snapPostD: after D.3 firm loan outflows applied
    vector<int64_t> snapPostABankLiq;
    vector<int64_t> snapPostCBankLiq;
    vector<int64_t> snapPostDBankLiq;
    for (int run = 0; run < runs; run++) {

        // Policy = 0 -> Bail-In
        // Policy = 1 -> Bail-Out
        for (int policy = 0; policy <= 1; policy++) {

            // Build the per-(run, policy) debug output folder once here so
            // all timestep files land in the same directory regardless of
            // how long the run takes. Rank 0 creates the folder; the path
            // is broadcast so every rank agrees on the name.
            std::string auditRunFolder;
            if (debug) {
                auditRunFolder = auditBuildRunFolder(
                    (unsigned int)mpiRank,
                    (unsigned int)mpiSize,
                    (unsigned int)run,
                    (unsigned int)policy
                );
            }

            // <!===========================================================!>
            // <!===================!Create Local Arrays!===================!>
            // <!===========================================================!>

            // Bank Arrays
            vector<int64_t> localBankLiquidity(bankCountForRank,
                initialBankLiquidity);
            vector<unsigned short> localBankInterestRate(bankCountForRank);
            vector<int64_t> wageDepositBuffer(bankCountTotal);
            vector<FirmLoanRequest> receivedLoanRequests(0);

            vector<vector<d2FirmLoanAcceptance>> firmLoanAcceptancesToRank(mpiSize);
            vector<uint64_t> bankFirmLoanOutflow(bankCountForRank);
            vector<d2FirmLoanAcceptance> receivedAcceptances;
            vector<vector<DebtEntry>> eLocalBankDebts(bankCountForRank);
            vector<unsigned int> bankDistressCount(bankCountForRank, 0);

            // Firm Arrays
            vector<int64_t> localFirmLiquidity(firmCountForRank,
                initialFirmLiquidity);
            vector<uint64_t> localFirmProductionCost(firmCountForRank,
                initialProductionCost);
            vector<vector<unsigned int>> localFirmNeighbors(firmCountForRank);
            vector<vector<DebtEntry>> localFirmDebts(firmCountForRank);
            vector<vector<FirmLoanRequest>> firmLoanRequestsToRank(mpiSize);

            // Worker Arrays
            vector<unsigned short> localWorkerBankID(workerCountForRank),
                localWorkerFirmID(workerCountForRank);

            // <!============================================================!>
            // <!==============!Initialize Worker Local Arrays!==============!>
            // <!============================================================!>

            // Randomly assign a bank to each worker
            // Randomly assign a firm to each worker
            for (int i = 0; i < localWorkerBankID.size(); i++) {
                localWorkerBankID[i] = selectForWorker(baseSeed, i +
                    workerGlobalStartIndex, bankCountTotal,
                    STREAM_WORKER_BANK_ASSIGN);

                localWorkerFirmID[i] = selectForWorker(baseSeed, i +
                    workerGlobalStartIndex, firmCountTotal,
                    STREAM_WORKER_FIRM_ASSIGN);
            }

            // <!============================================================!>
            // <!===============!Initialize Local Bank Arrays!===============!>
            // <!============================================================!>

            for (unsigned int i = 0; i < bankCountForRank; i++) {
                // Get the global bank ID
                unsigned int bankGlobalId = bankGlobalStartIndex + i;

                // Find the seed for this bank's interest rate
                uint64_t bankSeed = makeSeed(
                    baseSeed,
                    run,
                    0,
                    bankGlobalId,
                    STREAM_INTEREST_RATE_SELECTION
                );

                localBankInterestRate[i] = randomInRangeFromSeed(
                    bankSeed,
                    minInterestRate,
                    maxInterestRate
                );
            }

            // <!===========================================================!>
            // <!==================!Build Interbank Graph!==================!>
            // <!===========================================================!>

            // Find the degree from density
            unsigned int degree = (unsigned int)(((uint64_t)(bankCountTotal
                - 1) * interbankDensity) / 100);

            // Clamp for safety. Uncomment if graph is exploding too fast
            // degree = min(degree, (unsigned int)maxInterbankLenderSamplingK);

            // Neighbor list for each local bank
            vector<vector<unsigned int>>
                localBankNeighbors(bankCountForRank);

            for (unsigned int localBorrower = 0; localBorrower <
                bankCountForRank; localBorrower++) {

                unsigned int borrowerGlobalId = bankGlobalStartIndex
                    + localBorrower;

                // Deterministic per run, policy, borrowerGlobalId
                uint64_t state = makeSeed(
                    baseSeed,
                    run,
                    0,
                    borrowerGlobalId,
                    STREAM_BANK_GRAPH_BUILD
                );

                // Sample degree # of unique neighbors excluding self
                localBankNeighbors[localBorrower].clear();

                if (degree == 0) {
                    continue;
                }

                // Create a temporary vector to track which banks have already been selected
                // Sized to total number of banks
                vector<uint8_t> used(bankCountTotal, 0);

                // Prevent self use by marking this borrower as already used
                // Makes sure a bank cannot be its own neighbor
                used[borrowerGlobalId] = 1;

                // Continue sampling until reach the degree as defined above
                while (localBankNeighbors[localBorrower].size() < degree) {

                    // Advance deterministic state to get next value
                    state = splitmix64Hash(state);

                    // Map random value into valid bank ID range
                    unsigned int candidate = (unsigned int)(state %
                        bankCountTotal);

                    // If this candidate bank has not already been selected
                    if (!used[candidate]) {
                        // Mark it as used
                        used[candidate] = 1;
                        // Add the candidate as a neighbor
                        localBankNeighbors[localBorrower].push_back(candidate);
                    }
                }
            }

            // <!============================================================!>
            // <!=================!Build Firm -> Bank Graph!=================!>
            // <!============================================================!>

            // Build firm graph 
            for (unsigned int localFirm = 0; localFirm < firmCountForRank; localFirm++) {
                unsigned int firmGlobalId = firmGlobalStartIndex + localFirm;

                buildFirmNeighbors(
                    baseSeed,
                    (unsigned int)run,
                    firmGlobalId,
                    bankCountTotal,
                    (unsigned int)firmLenderDegree,
                    STREAM_FIRM_BANK_GRAPH_BUILD,
                    localFirmNeighbors[localFirm]
                );
            }

            // <!============================================================!>
            // <!================!Pre-compute Workforce Info!================!>
            // <!============================================================!>

            // How much each firm has to pay in total
            vector<uint64_t> localFirmWorkforceCost(firmCountForRank, 0);

            a0ComputeFirmWorkforceCost(
                firmCountTotal,
                mpiRank,
                mpiSize,
                firmGlobalStartIndex,
                firmCountForRank,
                localWorkerFirmID,
                wage,
                localFirmWorkforceCost
            );

            // <!============================================================!>
            // <!====================!Main TimeStep Loop!====================!>
            // <!============================================================!>

            for (int step = 0; step < timesteps; step++) {

                //..............PHASE A: Firms and Worker Updates..............

                // Firms compute loans before worker deposits are settled,
                // which may cause them to miss out on a bank they otherwise
                // could've used as a candidate, modeling short-term 
                // overreaction to liquidity shocks
                vector<int64_t> bankIncomingFromFirmRepay(bankCountTotal, 0);
                for (int f = 0; f < firmCountForRank; f++) {
                    // A.1: Find the random shock value
                    uint64_t shockSeed = makeSeed(
                        baseSeed,
                        run,
                        step,
                        f + firmGlobalStartIndex,
                        STREAM_SHOCK_GENERATION
                    );

                    short shock = randomInRangeFromSeed(
                        shockSeed,
                        shockMultiplierMin,
                        shockMultiplierMax
                    );

                    bool isNegative = (shockSeed & 1ULL);

                    if (isNegative) {
                        shock *= -1;
                    }

                    // A.2: Affect the production cost using the shock value

                    localFirmProductionCost[f] *= (1 + (shock / 100.0));

                    // A.3: Find the random profit multiplier and affect profit

                    uint64_t profitSeed = makeSeed(
                        baseSeed,
                        run,
                        step,
                        f + firmGlobalStartIndex,
                        STREAM_PROFIT_MULTIPLIER_GENERATION
                    );

                    short profitMultiplier = randomInRangeFromSeed(
                        profitSeed,
                        profitMultiplierMin,
                        profitMultiplierMax
                    );

                    // A.4: Apply profit and take away workforce cost

                    localFirmLiquidity[f] += localFirmProductionCost[f] *
                        (profitMultiplier / 100.0) - localFirmWorkforceCost[f];

                    // A.5: Repay bank loans (firmRepayPercent % of the loan)

                    // If liquidity is smaller than the percentage of the loan
                    // -> Pay all of liquidity

                    // If liquidity is negative
                    // -> Don't pay at all

                    repayFirmLoansProRata(
                        localFirmLiquidity[f],
                        localFirmDebts[f],
                        firmRepayPercent,
                        bankIncomingFromFirmRepay
                    );
                }

                // A.5.5: Bank pays employees their wage.
                // They keep a percentage of it in their bank accounts
                // A bank employee deposits into the bank they work for
                for (int bankID = 0; bankID < bankCountForRank; bankID++) {
                    localBankLiquidity[bankID] -=
                        (int64_t)bankEmployeeWage *
                        (int64_t)bankWorkerCount *
                        (int64_t)wageConsumptionPercent / 100;
                }

                // Snapshot A: bank liquidity after employee wages, before C deposits
                if (debug && step % logEvery == 0)
                    snapPostABankLiq = localBankLiquidity;

                // A.6: Request loans if needed

                // Clear request buffer
                for (int r = 0; r < mpiSize; r++) {
                    firmLoanRequestsToRank[r].clear();
                }

                a6RequestFirmLoans(
                    run,
                    baseSeed,
                    bankCountTotal,
                    mpiSize,
                    minInterestRate,
                    maxInterestRate,
                    firmGlobalStartIndex,
                    localFirmLiquidity,
                    localFirmNeighbors,
                    firmLoanRequestsToRank
                );

                //.....................PHASE B: Worker Deposits....................

                // B.1: Consume part of wage and deposit the rest
                // (record transaction to buffer)
                uint64_t depositPerW = ((uint64_t)wage *
                    (uint64_t)(100 - wageConsumptionPercent)) / 100ULL;

                // Go through all workers' banks and add their remaining wage
                // to the buffer

                for (int wID = 0; wID < localWorkerBankID.size(); wID++) {
                    wageDepositBuffer[localWorkerBankID[wID]] += depositPerW;
                }

                //...................PHASE C: Send Money to Banks .................

                // C.1: Send worker money and firm loan repayment deltas to the 
                // corresponding bank (account for global and local ID conversion)
                // C.2: Update bank balances
                // C.1 and C2 were bundled together

                c1SendMoneyToBanks(
                    bankCountTotal,
                    (unsigned int)mpiRank,
                    (unsigned int)mpiSize,
                    bankGlobalStartIndex,
                    bankCountForRank,
                    localBankLiquidity,
                    wageDepositBuffer,
                    bankIncomingFromFirmRepay
                );

                // Snapshot C: bank liquidity after worker deposits + firm repayments received
                if (debug && step % logEvery == 0)
                    snapPostCBankLiq = localBankLiquidity;

                //......................PHASE D: Request Loans.....................

                // Clear before starting the requesting process
                receivedLoanRequests.clear();

                // D.1: Send loan request buffers to the banks
                d1SendFirmLoanRequests(
                    mpiSize,
                    firmLoanRequestsToRank,
                    receivedLoanRequests
                );

                // D.2: Banks process loan requests by ID
                d2ProcessFirmLoanRequests(
                    firmCountTotal,
                    bankCountForRank,
                    bankGlobalStartIndex,
                    mpiSize,
                    maxFirmLoanPercent,
                    localBankLiquidity,
                    receivedLoanRequests,
                    firmLoanAcceptancesToRank,
                    bankFirmLoanOutflow
                );

                // D.3: Banks update their liquidity to reflect accepted loans
                for (unsigned int b = 0; b < bankCountForRank; b++) {
                    localBankLiquidity[b] -= (int64_t)bankFirmLoanOutflow[b];

                    // Liquidity should not go negative from lending
                    if (bankFirmLoanOutflow[b] > 0 && localBankLiquidity[b] < 0) {
                        printf("ERROR: Bank liquidity negative after D3\n");
                    }
                }

                // Snapshot D: bank liquidity after firm loan outflows (pre-E, pre-F)
                if (debug && step % logEvery == 0)
                    snapPostDBankLiq = localBankLiquidity;

                // D.4: Banks send back acceptances 
                d4SendFirmLoanAcceptances(
                    mpiSize,
                    firmLoanAcceptancesToRank,
                    receivedAcceptances
                );

                // D.5: If accepted, firms update liquidity and loan list
                d5ApplyFirmLoanAcceptances(
                    firmGlobalStartIndex,
                    localFirmLiquidity,
                    localFirmDebts,
                    receivedAcceptances
                );

                //...................PHASE E: Bank Updates..................

                // E.1: Repay interbank loans (bankRepayPercent % of them)
                // The repayments are saved in a buffer
                // The subtraction is applied instantly

                e1RepayInterbankLoans(
                    bankCountTotal,
                    mpiRank,
                    mpiSize,
                    bankGlobalStartIndex,
                    bankCountForRank,
                    bankRepayPercent,
                    localBankLiquidity,
                    eLocalBankDebts
                );

                // E.2: Accrue interest on outstanding loans
                // Applied by the borrower
                // Iterate debt lists for firms and banks, and update loan amount
                e2ApplyInterestOnAllLoans(
                    eLocalBankDebts,
                    localFirmDebts,
                    baseSeed,
                    run,
                    minInterestRate,
                    maxInterestRate
                );

                // E.3: If the bank's liquidity is negative, attempt to 
                // borrow from another bank (by sending a request message)
                e3BorrowInterbankIfNegativeLiquidity(
                    bankCountTotal,
                    mpiRank,
                    mpiSize,
                    bankGlobalStartIndex,
                    bankCountForRank,
                    maxInterbankLenderSamplingK,
                    maxInterbankLoanPercent,
                    baseSeed,
                    run,
                    minInterestRate,
                    maxInterestRate,
                    localBankNeighbors,
                    localBankLiquidity,
                    eLocalBankDebts
                );

                if (debug && step % logEvery == 0) {
                    preFBankLiquidity = localBankLiquidity;
                    preFBankDebts = eLocalBankDebts;
                    preFBankDistressCount = bankDistressCount;
                }

                //...................PHASE F: Insolvency + Policy..................

                // F.1: Insolvency check for each local bank
                // Insolvent means bank's liquidity is smaller than the total debt
                // If insolvent, increment the bank's distress counter
                // If solvent, reset the distress counter to 0
                // F.2: If insolvent and the bank has been distressed for at least
                // interventionDelay timesteps, apply policy (once per distress episode)
                //
                // Policy 0 (Bail-In):
                //  Deficit = amount needed
                //  L = total interbank liabilities (debts)
                //  haircutRatio = min(1.0, deficit / L)
                //  haircutFraction = haircutRatio * bailInCoveragePercent / 100.0
                //
                //
                // Policy 1 (Bail-Out, interbank only):
                //  Each interbank loan owed by this bank is reduced multiplicatively
                //  by (1 - haircutFraction).
                //  Deficit = amount needed to restore solvency
                //  Injection = Deficit * bailOutCoveragePercent / 100.0
                //
                //  The injection is added directly to the bank's liquidity.   
                //
                // After the policy is applied, mark the intervention as applied for
                // this distress episode so it is not applied again until the bank
                // returns to solvency.
                //
                // For both policies, track how much bail was granted this step
                // (cash injected for bail-out or total debt removed for bail-in)
                uint64_t bailMoneyThisStep = 0, graceMoneyThisStep = 0;
                f1f2UpdateBankDistressAndApplyIntervention(
                    bankCountForRank,
                    interventionDelay,
                    policy,
                    bailInCoveragePercent,
                    bailOutCoveragePercent,
                    localBankLiquidity,
                    eLocalBankDebts,
                    bankDistressCount,
                    bailMoneyThisStep,
                    graceMoneyThisStep
                );

                // F.3: Update total bail money
                // After the bail money for this timestep is finalized,
                // update the total tracker
                if (policy == 0) {
                    totalBailInMoney += bailMoneyThisStep;
                }
                else {
                    totalBailOutMoney += bailMoneyThisStep;
                }

                // F.6: Update total grace money
                // After all banks have been processed, add the grace money
                // used this timestep to the total grace money tracker
                totalGraceMoney = graceMoneyThisStep;

                if (debug && step % logEvery == 0) {
                    JSONL_AUDIT_TIMESTEP(
                        (unsigned int)run,
                        (unsigned int)step,
                        (unsigned int)policy,
                        (unsigned int)mpiRank,
                        (unsigned int)mpiSize,
                        auditRunFolder,
                        bankCountTotal,
                        bankGlobalStartIndex,
                        bankCountForRank,
                        firmCountTotal,
                        firmGlobalStartIndex,
                        firmCountForRank,
                        workerCountForRank,
                        workerGlobalStartIndex,
                        // Phase snapshots (A, C, D)
                        snapPostABankLiq,
                        snapPostCBankLiq,
                        snapPostDBankLiq,
                        // Pre-F snapshots
                        preFBankLiquidity,
                        preFBankDebts,
                        preFBankDistressCount,
                        // Post-F live state
                        localBankLiquidity,
                        eLocalBankDebts,
                        bankDistressCount,
                        // Phase F outputs
                        bailMoneyThisStep,
                        graceMoneyThisStep,
                        // Firm state
                        localFirmLiquidity,
                        localFirmDebts,
                        localBankNeighbors,
                        localFirmNeighbors,
                        localFirmProductionCost,
                        // Simulation parameters
                        baseSeed,
                        bankWorkerCount,
                        bankEmployeeWage,
                        wage,
                        workerCountTotal,
                        localWorkerBankID,
                        localWorkerFirmID,
                        minInterestRate,
                        maxInterestRate,
                        shockMultiplierMin,
                        shockMultiplierMax,
                        profitMultiplierMin,
                        profitMultiplierMax,
                        wageConsumptionPercent,
                        localFirmWorkforceCost
                    );
                }
            }
        }
    }

}