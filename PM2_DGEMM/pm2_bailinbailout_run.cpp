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

    /* System Parameters */
    const unsigned int initialBankLiquidity,
    const unsigned int initialFirmLiquidity,
    const unsigned int initialProductionCost,
    const unsigned int wage,

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
    unsigned short maxInterbankLoanPercent,
    unsigned short maxFirmLoanPercent,
    unsigned short firmLenderDegree,
    unsigned short firmRepayPercent,
    unsigned short bankRepayPercent,

    /* Bail-In Parameters */
    unsigned short bailInCoveragePercent,

    /* Bail-Out Parameters */
    unsigned short bailOutCoveragePercent
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
        bankRepayPercent, bailInCoveragePercent, bailOutCoveragePercent)) {
        return;
    }

    // <!===================================================================!>
    // <!=========================!MPI Information!=========================!>
    // <!===================================================================!>

    int mpiRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

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

    for (int run = 0; run < runs; run++) {
        // Policy = 0 -> Bail-In
        // Policy = 1 -> Bail-Out
        for (int policy = 0; policy <= 1; policy++) {

            // <!===========================================================!>
            // <!===================!Create Local Arrays!===================!>
            // <!===========================================================!>

            // Bank Arrays
            vector<int64_t> localBankLiquidity(bankCountForRank,
                initialBankLiquidity),
                localBankDebt(bankCountForRank, 0);
            vector<vector<unsigned int>> localBankDebtsOwedToBankID(
                bankCountForRank, vector<unsigned int>());
            vector<unsigned short> localBankInterestRate(bankCountForRank);
            vector<int64_t> wageDepositBuffer(bankCountTotal);
            vector<FirmLoanRequest> receivedLoanRequests(0);

            // Firm Arrays
            vector<int64_t> localFirmLiquidity(initialFirmLiquidity);
            vector<uint64_t> localFirmProductionCost(initialProductionCost);
            vector<vector<unsigned int>> localFirmNeighbors(firmCountForRank);
            vector<vector<DebtEntry>> localFirmDebts(firmCountForRank);
            vector<vector<FirmLoanRequest>> firmLoanRequestsToRank(mpiSize);

            // Worker Arrays
            vector<unsigned short> localWorkerBankID, localWorkerFirmID;

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

                // FInd the seed for this bank's interest rate
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
                std::vector<uint8_t> used(bankCountTotal, 0);

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

            for (unsigned short& worker : localWorkerFirmID) {
                localFirmWorkforceCost[localWorkerFirmID[worker]] += wage;
            }

            // <!============================================================!>
            // <!====================!Main TimeStep Loop!====================!>
            // <!============================================================!>

            for (int step = 0; step < timesteps; step++) {
                //......................PHASE A: Firm Updates......................

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

                // A.6: Request loans if needed

                // If the liquidity is negative, request a loan
                // Sample minimum of firmLenderSampleK and degree # of 
                // candidate banks from the adjacency list
                // If lender has enough liquidity and lower interest rate 
                // than the currently best choice, set as the best choice
                // Save the request in a buffer

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

                // D.3: Banks update their liquidity to reflect accepted loans

                // D.4: Banks send back acceptances 

                // D.5: If accepted, firms update liquidity and loan list

                //...................PHASE E: Bank Updates..................

                    // E.1: Repay interbank loans (bankRepayPercent % of them)
                    // The repayments are saved in a buffer
                    // The subtraction is applied instantly

                    // E.2: Accrue interest on outstanding loans
                    // Applied by the borrower
                    // Iterate debt lists for firms and banks, and update loan amount

                    // E.3: If the bank's liquidity is negative, attempt to 
                    // borrow from another bank (by sending a request message)
                    // Sample up to maxInterbankLenderSamplingK from the adjacency list

                    // Enforce maxInterbankLoanPercent capacity constraint
                    // Send loan requests to candidate lender bank owners
                    // Lender banks process requests deterministically and grant up to capacity
                    // If granted, write a delta for:
                    // Which global bank ID accepted their request
                    // How much this bank owes them

                //...................PHASE F: Settle Interbank Money.................

                    // F.1: Send interbank repayment deltas (and any loan disbursement deltas)
                    // to the corresponding bank owners

                    // F.2: Bank owners update their liquidity 

                    // F.3: Clear interbank buffers

                //...................PHASE G: Insolvency + Policy..................

                    // G.1: Insolvency check for each local bank
                    // Insolvent means bank's liquidity is smaller than the total debt

                    // G.2: If insolvent, apply policy
                    // Policy 0 (Bail-Out):
                    //   Deficit = amount needed to restore solvency
                    //   Injection = Deficit * bailOutCoveragePercent/100.0
                    //
                    // Policy 1 (Bail-In, interbank only):
                    //   Deficit = amount needed
                    //   L = total interbank liabilities (debts)
                    //   haircutRatio = min(1.0, deficit / L)
                    //   haircut fraction = haircutRatio * bailInCoveragePercent/100.0
                    //   Each interbank loan owed by this bank gets lowered by the
                    //   haircut fraction (multiplicatively)

                    // G.3: Track bail
                    // For both policies, track how much bail was granted this step

                    // G.4: Update total bail money
                    // After the bail money for this timestep is finalized,
                    // update the total tracker

                    // G.5: Grace money for failure
                    // If the coverege is not 100%, the bank might fail
                    // Find the remaining deficit after policy application
                    // Add the deficit to the grace money tracker for this timestep

                    // G.6: Update total grace money
                    // After all banks have been rescued, add the grace money
                    // used this time step to the total grace money tracker

            }
        }
    }

}