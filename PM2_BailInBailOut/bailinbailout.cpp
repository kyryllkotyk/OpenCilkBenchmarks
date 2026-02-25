#include "bailinbailout.h"

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

            // Firm Arrays
            vector<int64_t> localFirmLiquidity(initialFirmLiquidity);
            vector<uint64_t> localFirmProductionCost(initialProductionCost);
            vector<vector<unsigned int>> localFirmNeighbors(firmCountForRank);
            vector<vector<DebtEntry>> localFirmDebts(firmCountForRank);

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

                    // A.6: Request loans if needed

                    // If the liquidity is negative, request a loan

                    // Sample minimum of firmLenderSampleK and degree # of 
                    // candidate banks from the adjacency list

                    // If lender has enough liquidity and lower interest rate 
                    // than the currently best choice, set as the best choice

                    // Save the request in a buffer


                }

                //.....................PHASE B: Worker Deposits....................

                    // B.1: Consume part of wage
                    // B.2: Deposit the income (record transaction to buffer)
                        // Buffer could be an array with index = bank ID
                        // that holds all the money that it should receive from
                        // the workers

                //...................PHASE C: Send Money to Banks .................

                    // C.1: Send worker money and firm loan repayment deltas to the 
                    // corresponding bank (account for global and local ID conversion)

                    // C.2: Update bank balances

                //......................PHASE D: Request Loans.....................

                    // D.1: Send loan request buffers to the banks

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




//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!HELPERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bool BailInBailOut::errorDetectionAndClamping(
    const unsigned int runs,
    const unsigned int timesteps,
    const unsigned int logEvery,
    const uint64_t baseSeed,

    const unsigned int bankCountTotal,
    const unsigned int firmCountTotal,
    const unsigned int workerCountTotal,

    const unsigned int initialBankLiquidity,
    const unsigned int initialFirmLiquidity,
    const unsigned int initialProductionCost,
    const unsigned int wage,

    unsigned short& wageConsumptionPercent,
    unsigned short& profitMultiplierMin,
    unsigned short& profitMultiplierMax,
    unsigned short& shockMultiplierMin,
    unsigned short& shockMultiplierMax,
    unsigned short& minInterestRate,
    unsigned short& maxInterestRate,

    unsigned short& interbankDensity,
    const unsigned short maxInterbankLenderSamplingK,
    unsigned short& maxInterbankLoanPercent,
    unsigned short& maxFirmLoanPercent,
    unsigned short& firmLenderDegree,
    unsigned short& firmRepayPercent,
    unsigned short& bankRepayPercent,

    unsigned short& bailInCoveragePercent,
    unsigned short& bailOutCoveragePercent
) {
    (void)baseSeed;

    // Ensure all multipliers are between 0 and 100
    // Helper to do that
    auto clampPercent = [](unsigned short value) -> unsigned short {
        return (value > 100) ? 100 : value;
        };

    wageConsumptionPercent = clampPercent(wageConsumptionPercent);
    profitMultiplierMin = clampPercent(profitMultiplierMin);
    profitMultiplierMax = clampPercent(profitMultiplierMax);
    shockMultiplierMin = clampPercent(shockMultiplierMin);
    shockMultiplierMax = clampPercent(shockMultiplierMax);

    interbankDensity = clampPercent(interbankDensity);
    maxInterbankLoanPercent = clampPercent(maxInterbankLoanPercent);
    maxFirmLoanPercent = clampPercent(maxFirmLoanPercent);

    firmRepayPercent = clampPercent(firmRepayPercent);
    bankRepayPercent = clampPercent(bankRepayPercent);

    bailInCoveragePercent = clampPercent(bailInCoveragePercent);
    bailOutCoveragePercent = clampPercent(bailOutCoveragePercent);

    // Ensure that the min & max make sense (max >= min)
    if (profitMultiplierMin > profitMultiplierMax) {
        unsigned short tmp = profitMultiplierMin;
        profitMultiplierMin = profitMultiplierMax;
        profitMultiplierMax = tmp;
    }

    if (shockMultiplierMin > shockMultiplierMax) {
        unsigned short tmp = shockMultiplierMin;
        shockMultiplierMin = shockMultiplierMax;
        shockMultiplierMax = tmp;
    }

    if (minInterestRate > maxInterestRate) {
        unsigned short tmp = minInterestRate;
        minInterestRate = maxInterestRate;
        maxInterestRate = tmp;
    }

    // <!=========================!Error Detection!=========================>

    // Out of bounds errors
    if (runs < 1) {
        cerr << "runs variable out of bounds (Must be at least 1)";
        return false;
    }

    if (timesteps < 1) {
        cerr << "timesteps variable out of bounds (Must be at least 1)";
        return false;
    }

    if (logEvery < 1) {
        cerr << "logEvery variable out of bounds (Must be at least 1)";
        return false;
    }

    if (bankCountTotal < 1) {
        cerr << "bankCountTotal variable out of bounds (Must be at least 1)";
        return false;
    }

    if (firmCountTotal < 1) {
        cerr << "firmCountTotal variable out of bounds (Must be at least 1)";
        return false;
    }

    if (workerCountTotal < 1) {
        cerr << "workerCountTotal variable out of bounds (Must be at least 1)";
        return false;
    }

    if (initialBankLiquidity < 1) {
        cerr << "initialBankLiquidity variable out of bounds (Must be at least 1)";
        return false;
    }

    if (initialFirmLiquidity < 1) {
        cerr << "initialBankLiquidity variable out of bounds (Must be at least 1)";
        return false;
    }

    if (initialProductionCost < 1) {
        cerr << "initialProductionCost variable out of bounds (Must be at least 1)";
        return false;
    }

    if (wage < 1) {
        cerr << "wage variable out of bounds (Must be at least 1)";
        return false;
    }

    if (maxInterbankLenderSamplingK < 1) {
        cerr << "maxInterbankLenderSamplingK variable out of bounds (Must be at least 1)";
        return false;
    }

    if (firmLenderDegree < 1) {
        cerr << "firmLenderSampleK variable out of bounds (Must be at least 1)";
        return false;
    }

    if (logEvery > timesteps) {
        cerr << "logEvery variable out of bounds (Cannot exceed timesteps)\n"
            << "Set to zero if you don't want logs";
        return false;
    }

    // Memory protection
    if (bankCountTotal > 10000000) {
        cerr << "bankCountTotal too large (Unreasonable allocation size)";
        return false;
    }

    if (timesteps > 1000000) {
        cerr << "timesteps too large (Not in scope of the simulation)";
        return false;
    }

    // Degenerate Cases
    if (wageConsumptionPercent == 0 &&
        interbankDensity == 0 &&
        maxInterbankLoanPercent == 0 &&
        maxFirmLoanPercent == 0 &&
        firmRepayPercent == 0 &&
        bankRepayPercent == 0 &&
        bailInCoveragePercent == 0 &&
        bailOutCoveragePercent == 0) {
        cerr << "configuration out of bounds (No transfers possible. Simulation will be degenerate)";
        return false;
    }

    return true;
}

pair<unsigned int, unsigned int> BailInBailOut::computeRange(
    unsigned int totalCount,
    unsigned int rank,
    unsigned int totalRanks) {

    unsigned int base = totalCount / totalRanks;
    unsigned int remainder = totalCount % totalRanks;

    pair<unsigned int, unsigned int> range;

    if (rank < remainder) {
        // Count
        range.first = base + 1;
        // Start
        range.second = rank * (base + 1);
    }
    else {
        range.first = base;
        range.second = remainder * (base + 1) + (rank - remainder) * base;
    }

    return range;
}

uint64_t BailInBailOut::Xoshiro256::splitmix64(uint64_t& x) {
    x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

BailInBailOut::Xoshiro256::Xoshiro256(uint64_t seed) {
    for (int i = 0; i < 4; i++) {
        s[i] = splitmix64(seed);
    }
}

uint64_t BailInBailOut::Xoshiro256::rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t BailInBailOut::Xoshiro256::next() {
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = rotl(s[3], 45);

    return result;
}

uint64_t BailInBailOut::Xoshiro256::nextInRange(uint64_t min, uint64_t max) {
    uint64_t range = max - min + 1;
    uint64_t x;
    uint64_t limit = std::numeric_limits<uint64_t>::max()
        - (std::numeric_limits<uint64_t>::max() % range);

    do {
        x = next();
    } while (x >= limit);

    return min + (x % range);
}

double BailInBailOut::Xoshiro256::nextDouble01() {
    return (next() >> 11) * (1.0 / 9007199254740992.0); // 2^53
}

double BailInBailOut::Xoshiro256::nextInRangeDouble(double min, double max) {
    return min + (max - min) * nextDouble01();
}


inline uint64_t BailInBailOut::splitmix64Hash(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

unsigned int BailInBailOut::selectForWorker(uint64_t baseSeed,
    unsigned int workerGlobalId,
    unsigned int countTotal,
    uint64_t streamTag
) {

    uint64_t seed = makeSeed(baseSeed, 0, 0, workerGlobalId, streamTag);
    BailInBailOut::Xoshiro256 rng(seed);

    // Bank IDs are global
    return (unsigned int)rng.nextInRange(0, (uint64_t)countTotal - 1);
}

inline uint64_t BailInBailOut::makeSeed(
    uint64_t baseSeed,
    unsigned int run,
    unsigned int timestep,
    unsigned int globalId,
    uint64_t streamTag
) {

    uint64_t key = baseSeed;

    // Mix in run number
    key ^= 0xD6E8FEB86659FD93ULL * (uint64_t)run;

    // Mix in timestep (0 for static parameters like graph/structure)
    key ^= 0xBF58476D1CE4E5B9ULL * (uint64_t)timestep;

    // Mix in entity global ID (bank/worker/firm)
    key ^= 0x9E3779B97F4A7C15ULL * (uint64_t)globalId;

    // Mix in stream tag to separate RNG 
    key ^= 0x94D049BB133111EBULL * streamTag;

    return splitmix64Hash(key);
}

inline unsigned short BailInBailOut::percent0to100FromMixedSeed(uint64_t seed) {
    // Assumes seed is already mixed
    uint64_t x = seed;

    const uint64_t range = 101;
    const uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

    // Retry by adding an odd constant 
    while (x >= limit) {
        x += 0x9e3779b97f4a7c15ULL;
    }

    return (unsigned short)(x % range);
}

inline uint64_t BailInBailOut::randomInRangeFromSeed(uint64_t seed, uint64_t minVal, uint64_t maxVal) {
    // Assumes mixed seed

    if (minVal > maxVal) {
        std::swap(minVal, maxVal);
    }

    uint64_t range = maxVal - minVal + 1;
    uint64_t x = seed;

    uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

    // Keep going until retry
    while (x >= limit) {
        x += 0x9e3779b97f4a7c15ULL;
    }

    return minVal + (x % range);
}

inline uint64_t BailInBailOut::percentFloorU64(uint64_t x, unsigned short pct) {
    return (x * (uint64_t)pct) / 100ULL;
}

// Pro-rata repayment for one firm.
// Writes lender receipts into bankIncomingFromFirmRepay for its global index
inline void BailInBailOut::repayFirmLoansProRata(
    int64_t& firmLiquidity,
    vector<DebtEntry>& debts,
    unsigned short firmRepayPercent,
    vector<int64_t>& bankIncomingFromFirmRepay
) {
    // Return if not enough liquidity, no debts, or simulation set to not repay
    if (firmLiquidity <= 0 || debts.empty() || firmRepayPercent == 0) {
        return;
    }

    // How many debts are owed
    const size_t n = debts.size();
    vector<uint64_t> scheduled(n, 0);

    uint64_t scheduledTotal = 0;

    // Go through all debts
    for (size_t i = 0; i < n; i++) {
        // How much this specific debt is
        uint64_t debtAmount = debts[i].amount;
        // If this debt is zero, return
        if (debtAmount == 0) {
            continue;
        }

        uint64_t s = percentFloorU64(debtAmount, firmRepayPercent);
        scheduled[i] = s;
        scheduledTotal += s;
    }

    // Nothing to pay
    if (scheduledTotal == 0) {
        return;
    }

    // Maximum money the firm can use to pay
    uint64_t canPay = (uint64_t)firmLiquidity;
    // How much they are allowed to pay total
    uint64_t payTotal = (scheduledTotal <= canPay) ? scheduledTotal : canPay;

    if (payTotal == 0) {
        return;
    }


    vector<uint64_t> paymentReceived(n, 0);
    vector<uint64_t> proportionalRemainder(n, 0);

    // Tracker of how much got assigned so far
    uint64_t sumPay = 0;
    for (size_t i = 0; i < n; i++) {
        if (scheduled[i] == 0) {
            continue;
        }

        // Payment is the floor of (payTotal * scheduled / scheduledTotal)
        unsigned __int128 scaledNumerator = (unsigned __int128)payTotal
            * (unsigned __int128)scheduled[i];
        uint64_t paymentAmount = (uint64_t)(scaledNumerator / scheduledTotal);
        uint64_t proportionalRemain = (uint64_t)(scaledNumerator %
            scheduledTotal);

        // Safety cap
        if (paymentAmount > debts[i].amount) {
            paymentAmount = debts[i].amount;
        }

        paymentReceived[i] = paymentAmount;
        proportionalRemainder[i] = proportionalRemain;
        sumPay += paymentAmount;
    }

    // Distribute leftover 1 by 1 by largest remainder 
    uint64_t leftover = payTotal - sumPay;
    if (leftover > 0) {
        vector<size_t> order;
        order.reserve(n);

        for (size_t i = 0; i < n; i++) {
            if (scheduled[i] == 0) {
                continue;
            }
            if (debts[i].amount > paymentReceived[i]) {
                order.push_back(i);
            }
        }

        sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            if (proportionalRemainder[a] != proportionalRemainder[b]) {
                return proportionalRemainder[a] > proportionalRemainder[b];
            }
            if (debts[a].lenderBankGlobalId != debts[b].lenderBankGlobalId) {
                return debts[a].lenderBankGlobalId < debts[b].lenderBankGlobalId;
            }
            return a < b;
        });

        for (size_t k = 0; k < order.size() && leftover > 0; k++) {
            size_t i = order[k];
            if (debts[i].amount > paymentReceived[i]) {
                paymentReceived[i] += 1;
                leftover -= 1;
            }
        }
    }

    // Apply payments
    uint64_t actuallyPaid = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t p = paymentReceived[i];
        if (p == 0) {
            continue;
        }

        debts[i].amount -= p;
        bankIncomingFromFirmRepay[debts[i].lenderBankGlobalId] += (int64_t)p;
        actuallyPaid += p;
    }

    firmLiquidity -= (int64_t)actuallyPaid;

    // Cleanup paid debts
    for (size_t i = 0; i < debts.size(); ) {
        if (debts[i].amount == 0) {
            debts[i] = debts.back();
            debts.pop_back();
        }
        else {
            i++;
        }
    }
}

inline void BailInBailOut::buildFirmNeighbors(
    uint64_t baseSeed,
    unsigned int run,
    unsigned int firmGlobalId,
    unsigned int bankCountTotal,
    unsigned int firmLenderDegree,
    uint64_t stream,
    vector<unsigned int>& outNeighbors
) {
    outNeighbors.clear();

    if (bankCountTotal == 0 || firmLenderDegree == 0) {
        return;
    }

    // Degree can't go above bankCountTotal
    unsigned int degree = std::min(firmLenderDegree, bankCountTotal);

    // Get the seed
    uint64_t seedy = makeSeed(
        baseSeed,
        run,
        0,
        firmGlobalId,
        stream
    );

    outNeighbors.reserve(degree);

    while (outNeighbors.size() < degree) {
        seedy = splitmix64Hash(seedy);

        unsigned int candidate = (unsigned int)(seedy % bankCountTotal);

        // Enforce uniqueness
        bool already = false;
        for (unsigned int x : outNeighbors) {
            if (x == candidate) {
                already = true;
                break;
            }
        }
        if (!already) {
            outNeighbors.push_back(candidate);
        }
    }
}

inline unsigned int BailInBailOut::bankOwnerRankFromGlobalId(
    unsigned int bankCountTotal,
    unsigned int totalRanks,
    unsigned int bankGlobalId
) {

    // Inverse mapping (global ID to owning rank, matches computeRange splitting)


    unsigned int base = bankCountTotal / totalRanks;
    unsigned int remainder = bankCountTotal % totalRanks;

    // Give +1 remainder to the first ranks until no remainder is left

    unsigned int cutoff = (base + 1) * remainder; 

    if (bankGlobalId < cutoff) {
        return bankGlobalId / (base + 1);
    }
    else {
        return remainder + (bankGlobalId - cutoff) / base;
    }
}

