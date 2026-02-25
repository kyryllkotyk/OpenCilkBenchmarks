#ifndef BAIL_IN_BAIL_OUT_
#define BAIL_IN_BAIL_OUT_ 


#include <iostream>
#include <vector>
using namespace std;

#define STREAM_GRAPH_BUILD 1ULL
#define STREAM_WORKER_BANK_ASSIGN  2ULL
#define STREAM_WORKER_FIRM_ASSIGN  3ULL
#define STREAM_INTEREST_RATE_SELECTION 4ULL
#define STREAM_SHOCK_GENERATION 5ULL
#define STREAM_PROFIT_MULTIPLIER_GENERATION 6ULL
class BailInBailOut
{
public:
    void runBenchmarkInAndOut(
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
        // Maximum percentage of the bank's liquidity it can offer as a loan
        unsigned short maxInterbankLoanPercent,
        unsigned short maxFirmLoanPercent,
        unsigned short firmLenderSampleK,
        unsigned short firmRepayPercent,
        unsigned short bankRepayPercent,

        /* Bail-In Parameters */
        unsigned short bailInCoveragePercent,

        /* Bail-Out Parameters */
        unsigned short bailOutCoveragePercent
    );

private:
    struct DebtEntry {
        // Global bank ID
        unsigned int lenderBankGlobalId; 
        // How much is owed
        uint64_t amount;
    };
    // Clamps fractions and checks for faulty parameters
    bool errorDetectionAndClamping(
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
        unsigned short& firmLenderSampleK,
        unsigned short& firmRepayPercent,
        unsigned short& bankRepayPercent,

        unsigned short& bailInCoveragePercent,
        unsigned short& bailOutCoveragePercent
    );

    // Finds range and starting index for the given rank
    pair<unsigned int, unsigned int> computeRange(
        unsigned int totalCount,
        unsigned int rank,
        unsigned int totalRanks
    );

    struct Xoshiro256 {
        uint64_t s[4];

        explicit Xoshiro256(uint64_t seed);

        static uint64_t splitmix64(uint64_t& x);
        static uint64_t rotl(uint64_t x, int k);

        uint64_t next();
        uint64_t nextInRange(uint64_t min, uint64_t max);

        double nextDouble01();
        double nextInRangeDouble(double min, double max);
    };

    static inline uint64_t splitmix64Hash(uint64_t x);

    unsigned int selectForWorker(
        uint64_t baseSeed,
        unsigned int workerGlobalId,
        unsigned int countTotal,
        uint64_t streamTag = 1
    );

    static inline uint64_t makeSeed(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int timestep,
        unsigned int globalId,
        uint64_t streamTag
    );

    uint64_t randomInRangeFromSeed(
        uint64_t seed,
        uint64_t minVal, 
        uint64_t maxVal);

    uint64_t percentFloorU64(
        uint64_t x,
        unsigned short pct
    );

    void repayFirmLoansProRata(
        int64_t& firmLiquidity, 
        vector<DebtEntry>& debts, 
        unsigned short firmRepayPercent,
        vector<int64_t>& bankIncomingFromFirmRepay
    );
};

#endif
