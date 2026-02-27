#ifndef BAIL_IN_BAIL_OUT_
#define BAIL_IN_BAIL_OUT_ 

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <limits>

#include <mpi.h>

using namespace std;

#define STREAM_BANK_GRAPH_BUILD 1ULL
#define STREAM_WORKER_BANK_ASSIGN  2ULL
#define STREAM_WORKER_FIRM_ASSIGN  3ULL
#define STREAM_INTEREST_RATE_SELECTION 4ULL
#define STREAM_SHOCK_GENERATION 5ULL
#define STREAM_PROFIT_MULTIPLIER_GENERATION 6ULL
#define STREAM_FIRM_BANK_GRAPH_BUILD 7ULL
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
        unsigned short firmLenderDegree,
        unsigned short firmRepayPercent,
        unsigned short bankRepayPercent,

        /* Bail-In Parameters */
        unsigned short bailInCoveragePercent,

        /* Bail-Out Parameters */
        unsigned short bailOutCoveragePercent,

        bool debug = false
    );

private:
    //<-----!!!!!!!!!!!!!!!!!!!!!!!!!STRUCTURES!!!!!!!!!!!!!!!!!!!!!!!!!----->

    struct DebtEntry {
        // Global bank ID
        unsigned int lenderBankGlobalId;
        // How much is owed
        uint64_t amount;
    };

    struct a0FirmWorkforceDeltaMessage {
        unsigned int firmGlobalId;
        uint64_t delta;
    };

    struct FirmLoanRequest {
        unsigned int firmGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountRequested;
        unsigned short lenderInterestRate;
    };

    struct c1BankDeltaMessage {
        unsigned int bankGlobalId;
        int64_t delta;
    };

    struct d2FirmLoanAcceptance {
        // Must use GLOBAL IDs
        unsigned int firmGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountGranted;
    };

    struct DDebugBankStats {
        unsigned int bankGlobalId;
        int64_t liquidityBeforeD3;
        uint64_t capThisStep;
        uint64_t outflowThisStep;
        int64_t liquidityAfterD3;

        uint64_t requestCount;
        uint64_t grantCount;
        uint64_t totalGranted;
    };

    struct DDebugFirmStats {
        unsigned int firmGlobalId;
        int64_t liquidityBeforeD5;
        int64_t liquidityAfterD5;
        uint64_t newDebtCount;
        uint64_t totalBorrowed;
    };

    struct DDebugGrantEdge {
        unsigned int firmGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountGranted;
    };


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

    struct e1InterbankRepaymentMessage {
        unsigned int lenderBankGlobalId;
        uint64_t amount;
    };

    //<-----!!!!!!!!!!!!!!!!!!!!!!!!MPI DATATYPES!!!!!!!!!!!!!!!!!!!!!!!!----->
    MPI_Datatype a0MakeFirmWorkforceDeltaType();
    MPI_Datatype d1MakeFirmLoanRequestType();
    MPI_Datatype d4MakeFirmLoanAcceptanceType();
    MPI_Datatype c1MakeBankDeltaMessageType();
    MPI_Datatype e1MakeInterbankRepaymentType();

    //TODO:: ORGANIZE BY FILE & PHASE
    // 
    // Ensures correctness up to D
    void ASSERT_CORRECTNESS_D(
        unsigned int bankCountTotal,
        unsigned short maxFirmLoanPercent,
        unsigned int bankCountForRank,
        vector<int64_t>& bankLiquidityBeforeD3,
        vector<uint64_t>& bankFirmLoanOutflow,
        vector<int64_t>& localBankLiquidityAfterD3,
        unsigned int firmCountForRank,
        vector<int64_t>& firmLiquidityBeforeD5,
        vector<int64_t>& localFirmLiquidityAfterD5,
        vector<unsigned int>& firmDebtCountBeforeD5,
        vector<vector<DebtEntry>>& localFirmDebtsAfterD5);

    void ASSERT_AND_LOG_D(
        unsigned int mpiRank, 
        unsigned int mpiSize,
        unsigned int bankCountTotal, 
        unsigned int bankGlobalStartIndex, 
        unsigned int bankCountForRank,
        unsigned int firmGlobalStartIndex, 
        unsigned int firmCountForRank,
        unsigned short maxFirmLoanPercent, 
        vector<int64_t>& bankLiquidityBeforeD3, 
        vector<int64_t>& firmLiquidityBeforeD5, 
        vector<unsigned int>& firmDebtCountBeforeD5, 
        vector<int64_t>& localBankLiquidity,
        vector<int64_t>& localFirmLiquidity, 
        vector<vector<DebtEntry>>& localFirmDebts,
        vector<uint64_t>& bankFirmLoanOutflow, 
        vector<DDebugGrantEdge>& grantsThisStep,
        bool printEdges
    );

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
        unsigned short& firmLenderDegree,
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


    static uint64_t splitmix64Hash(
        uint64_t x
    );

    unsigned int selectForWorker(
        uint64_t baseSeed,
        unsigned int workerGlobalId,
        unsigned int countTotal,
        uint64_t streamTag = 1
    );

    static uint64_t makeSeed(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int timestep,
        unsigned int globalId,
        uint64_t streamTag
    );

    unsigned short percent0to100FromMixedSeed(
        uint64_t seed
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

    void buildFirmNeighbors(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int firmGlobalId,
        unsigned int bankCountTotal,
        unsigned int firmLenderDegree,
        uint64_t stream,
        vector<unsigned int>& outNeighbors);

    unsigned int ownerRankFromGlobalId(
        unsigned int bankCountTotal,
        unsigned int totalRanks,
        unsigned int globalId
    );

    unsigned int getInterestRate(
        uint64_t baseSeed, 
        unsigned int run,
        unsigned int globalId, 
        uint64_t minVal,
        uint64_t maxVal
    );

    void a0ComputeFirmWorkforceCost(
        unsigned int firmCountTotal,
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int firmGlobalStartIndex,
        unsigned int firmCountForRank,
        const vector<unsigned short>& localWorkerFirmID,
        unsigned int wage,
        vector<uint64_t>& localFirmWorkforceCost
    );

    void a6RequestFirmLoans(
        unsigned int run,
        uint64_t baseSeed,
        unsigned int bankCountTotal,
        unsigned int mpiSize,
        unsigned short minInterestRate,
        unsigned short maxInterestRate,
        unsigned int firmGlobalStartIndex,
        vector<int64_t>& localFirmLiquidity,
        vector<vector<unsigned int>>& localFirmNeighbors,
        vector<vector<FirmLoanRequest>>& firmLoanRequestsToRank
    );

    void c1SendMoneyToBanks(
        unsigned int bankCountTotal,
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int bankGlobalStartIndex,
        unsigned int bankCountForRank,
        vector<int64_t>& localBankLiquidity,
        vector<int64_t>& wageDepositBuffer,
        vector<int64_t>& bankIncomingFromFirmRepay
    );

    void d1SendFirmLoanRequests(
        unsigned int mpiSize,
        vector<vector<FirmLoanRequest>>& firmLoanRequestsToRank,
        vector<FirmLoanRequest>& receivedLoanRequests
    );

    void d2ProcessFirmLoanRequests(
        unsigned int firmCountTotal,
        unsigned int bankCountForRank,
        unsigned int bankGlobalStartIndex,
        unsigned int mpiSize,
        unsigned short maxFirmLoanPercent,
        vector<int64_t>& localBankLiquidity,
        vector<FirmLoanRequest>& receivedLoanRequests,
        vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank,
        vector<uint64_t>& bankFirmLoanOutflow,

        vector<DDebugGrantEdge>* grantsThisStep,
        bool debug = false
    );

    void d4SendFirmLoanAcceptances(
        unsigned int mpiSize, 
        vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank, 
        vector<d2FirmLoanAcceptance>& receivedAcceptances
    );

    void d5ApplyFirmLoanAcceptances(
        unsigned int firmGlobalStartIndex,
        vector<int64_t>& localFirmLiquidity, 
        vector<vector<DebtEntry>>& localFirmDebts, 
        vector<d2FirmLoanAcceptance>& receivedAcceptances
    );

    void e1RepayInterbankLoans(
        unsigned int bankCountTotal, 
        unsigned int mpiRank,
        unsigned int mpiSize, 
        unsigned int bankGlobalStartIndex, 
        unsigned int bankCountForRank,
        unsigned short bankRepayPercent,
        vector<int64_t>& localBankLiquidity,
        vector<vector<DebtEntry>>& localBankDebts
    );


    // Exchange all helper
    template<typename T>
    void allToAllvExchange(
        unsigned int mpiSize,
        vector<vector<T>>& toRank,
        vector<T>& received,
        MPI_Datatype mpiType
    ) {
        vector<int> sendCounts(mpiSize, 0);
        for (unsigned int r = 0; r < mpiSize; r++) {
            sendCounts[r] = (int)toRank[r].size();
        }

        vector<int> recvCounts(mpiSize, 0);

        MPI_Alltoall(
            sendCounts.data(), 1, MPI_INT,
            recvCounts.data(), 1, MPI_INT,
            MPI_COMM_WORLD
        );

        vector<int> sendDispls(mpiSize, 0);
        vector<int> recvDispls(mpiSize, 0);

        int sendTotal = 0;
        int recvTotal = 0;

        for (unsigned int r = 0; r < mpiSize; r++) {
            sendDispls[r] = sendTotal;
            sendTotal += sendCounts[r];

            recvDispls[r] = recvTotal;
            recvTotal += recvCounts[r];
        }

        vector<T> sendBuf;
        sendBuf.reserve((unsigned int)sendTotal);

        for (unsigned int r = 0; r < mpiSize; r++) {
            for (unsigned int i = 0; i < toRank[r].size(); i++) {
                sendBuf.push_back(toRank[r][i]);
            }
        }

        received.clear();
        received.resize((unsigned int)recvTotal);

        MPI_Alltoallv(
            sendBuf.data(),
            sendCounts.data(),
            sendDispls.data(),
            mpiType,
            received.data(),
            recvCounts.data(),
            recvDispls.data(),
            mpiType,
            MPI_COMM_WORLD
        );

        // Clear outboxes
        for (unsigned int r = 0; r < mpiSize; r++) {
            toRank[r].clear();
        }
    }
};

#endif