#ifndef BAIL_IN_BAIL_OUT_
#define BAIL_IN_BAIL_OUT_ 

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <limits>
#include <cstdio> // For debug

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

    struct e3InterbankLoanRequest {
        unsigned int borrowerBankGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountRequested;
        unsigned short lenderInterestRate;
    };

    struct e3InterbankLoanAcceptance {
        unsigned int borrowerBankGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountGranted;
    };
    
    struct e3NeighborLiquidityMessage {
        unsigned int bankGlobalId;
        int64_t liquidity;
    };

    struct e3LenderOption {
        unsigned int lenderGlobalId;
        unsigned short interestRate;
        int64_t liquidity;
    };

    //<-----!!!!!!!!!!!!!!!!!!!!!!!!MPI DATATYPES!!!!!!!!!!!!!!!!!!!!!!!!----->
    MPI_Datatype a0MakeFirmWorkforceDeltaType();
    MPI_Datatype d1MakeFirmLoanRequestType();
    MPI_Datatype d4MakeFirmLoanAcceptanceType();
    MPI_Datatype c1MakeBankDeltaMessageType();
    MPI_Datatype e1MakeInterbankRepaymentType();
    MPI_Datatype e3MakeInterbankLoanRequestType();
    MPI_Datatype e3MakeInterbankLoanAcceptanceType();
    MPI_Datatype e3MakeNeighborLiquidityMessageType();
    
    //<-----!!!!!!!!!!!!!!!!!!!!!!!!DEBUG PROGRAM!!!!!!!!!!!!!!!!!!!!!!!!----->
    struct PreFVerifySnapshots {
        // Before E.1
        vector<int64_t> bankLiquidityBeforeE1;
        vector<vector<DebtEntry>> bankDebtsBeforeE1;
        vector<int64_t> firmLiquidityBeforeE1;
        vector<vector<DebtEntry>> firmDebtsBeforeE1;

        // After E.1 (before E.2)
        vector<int64_t> bankLiquidityAfterE1;
        vector<vector<DebtEntry>> bankDebtsAfterE1;
        vector<int64_t> firmLiquidityAfterE1;
        vector<vector<DebtEntry>> firmDebtsAfterE1;

        // After E.2 (before E.3)
        vector<int64_t> bankLiquidityAfterE2;
        vector<vector<DebtEntry>> bankDebtsAfterE2;
        vector<int64_t> firmLiquidityAfterE2;
        vector<vector<DebtEntry>> firmDebtsAfterE2;

        // After E.3 (before Phase F)
        vector<int64_t> bankLiquidityAfterE3;
        vector<vector<DebtEntry>> bankDebtsAfterE3;
        vector<int64_t> firmLiquidityAfterE3;
        vector<vector<DebtEntry>> firmDebtsAfterE3;
    };

    static void gatherAndPrintReportRank0Only(
        unsigned int mpiRank,
        unsigned int mpiSize,
        const std::string& localReport
    );

    static uint64_t sumDebtU64(
        const vector<BailInBailOut::DebtEntry>& debts
    );

    void ASSERT_AND_REPORT_PRE_F(
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int run,
        unsigned int timestep,
        unsigned int bankCountTotal,
        unsigned int bankGlobalStartIndex,
        unsigned int bankCountForRank,
        unsigned int firmCountTotal,
        unsigned int firmGlobalStartIndex,
        unsigned int firmCountForRank,
        vector<int64_t>& localBankLiquidity,
        vector<vector<DebtEntry>>& localBankDebts,
        vector<int64_t>& localFirmLiquidity,
        vector<vector<DebtEntry>>& localFirmDebts,
        vector<vector<unsigned int>>& localBankNeighbors,
        uint64_t baseSeed,
        unsigned int bankWorkerCount,
        unsigned int bankEmployeeWage,
        unsigned int wage,
        unsigned int workerCountTotal,
        unsigned short minInterestRate,
        unsigned short maxInterestRate,
        unsigned short shockMultiplierMin,
        unsigned short shockMultiplierMax,
        unsigned short profitMultiplierMin,
        unsigned short profitMultiplierMax,
        unsigned short wageConsumptionPercent,
        bool printAllBanks,
        bool hardAbortOnFailure
    );

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

    void e2ApplyInterestOnAllLoans(
        vector<vector<DebtEntry>>& eLocalBankDebts,
        vector<vector<DebtEntry>>& localFirmDebts,
        uint64_t baseSeed,
        unsigned int run,
        uint64_t minVal,
        uint64_t maxVal
    );

    void e3BorrowInterbankIfNegativeLiquidity(
        unsigned int bankCountTotal,
        unsigned int mpiRank,
        unsigned int mpiSize,
        unsigned int bankGlobalStartIndex,
        unsigned int bankCountForRank,
        unsigned short maxInterbankLenderSamplingK,
        unsigned short maxInterbankLoanPercent,
        uint64_t baseSeed,
        unsigned int run,
        unsigned short minInterestRate,
        unsigned short maxInterestRate,
        vector<vector<unsigned int>>& localBankNeighbors,
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

    template<typename T>
    void sparseDirectExchange(
        unsigned int mpiRank,
        unsigned int mpiSize,
        vector<vector<T>>& toRank,
        vector<T>& received,
        MPI_Datatype mpiType,
        int countTag,
        int dataTag
    ) {

        // Indicates whether THIS rank will send any payload to dst rank
        vector<int> willSendToRank(mpiSize, 0);
        for (unsigned int dstRank = 0; dstRank < mpiSize; dstRank++) {
            willSendToRank[dstRank] = (toRank[dstRank].empty() ? 0 : 1);
        }

        // Build a global send-mask matrix: [srcRank][dstRank]
        vector<int> allSendMasks(mpiSize * mpiSize, 0);

        MPI_Allgather(
            willSendToRank.data(),
            (int)mpiSize,
            MPI_INT,
            allSendMasks.data(),
            (int)mpiSize,
            MPI_INT,
            MPI_COMM_WORLD
        );

        // Indicates whether THIS rank will receive any payload from src rank
        vector<int> willReceiveFromRank(mpiSize, 0);
        for (unsigned int srcRank = 0; srcRank < mpiSize; srcRank++) {
            willReceiveFromRank[srcRank] =
                allSendMasks[(size_t)srcRank * (size_t)mpiSize + (size_t)mpiRank];
        }

        // -----------------------------
        // Exchange counts
        // -----------------------------

        vector<int> incomingCounts(mpiSize, 0);
        vector<MPI_Request> countReceives;

        for (unsigned int srcRank = 0; srcRank < mpiSize; srcRank++) {
            if (srcRank == mpiRank) {
                continue;
            }
            if (!willReceiveFromRank[srcRank]) {
                continue;
            }

            /*
            int count = (int)toRank[srcRank].size();
            if (count <= 0) {
                continue;
            }
            size_t elemSize = sizeof(T);
            size_t totalBytes = elemSize * (size_t)count;
            std::printf("[Rank %u] RECEIVE to rank %u: count=%d elemSize=%zu totalBytes=%zu countTag=%d dataTag=%d\n",
                mpiRank, srcRank, count, elemSize, totalBytes, countTag, dataTag);
            std::fflush(stdout);
            */

            MPI_Request req;
            MPI_Irecv(
                &incomingCounts[srcRank],
                1,
                MPI_INT,
                (int)srcRank,
                countTag,
                MPI_COMM_WORLD,
                &req
            );
            countReceives.push_back(req);
        }

        vector<MPI_Request> countSends;

        for (unsigned int dstRank = 0; dstRank < mpiSize; dstRank++) {
            if (dstRank == mpiRank) {
                continue;
            }
            if (!willSendToRank[dstRank]) {
                continue;
            }

            int outCount = (int)toRank[dstRank].size();

            /*
            int count = (int)toRank[dstRank].size();
            if (count <= 0) {
                continue;
            }

            size_t elemSize = sizeof(T);
            size_t totalBytes = elemSize * (size_t)count;
            printf("[Rank %u] SEND to rank %u: count=%d elemSize=%zu totalBytes=%zu countTag=%d dataTag=%d\n",
                (unsigned)mpiRank,
                (unsigned)dstRank,
                (int)count,
                (size_t)elemSize,
                (size_t)totalBytes,
                (int)countTag,
                (int)dataTag);
            fflush(stdout);
            */

            MPI_Request req;
            MPI_Isend(
                &outCount,
                1,
                MPI_INT,
                (int)dstRank,
                countTag,
                MPI_COMM_WORLD,
                &req
            );
            countSends.push_back(req);
        }

        if (!countReceives.empty()) {
            MPI_Waitall((int)countReceives.size(), countReceives.data(), MPI_STATUSES_IGNORE);
        }
        if (!countSends.empty()) {
            MPI_Waitall((int)countSends.size(), countSends.data(), MPI_STATUSES_IGNORE);
        }

        // -----------------------------
        // Exchange payloads
        // -----------------------------

        vector<vector<T>> incomingByRank(mpiSize);
        vector<MPI_Request> payloadReceives;

        for (unsigned int srcRank = 0; srcRank < mpiSize; srcRank++) {
            if (srcRank == mpiRank) {
                continue;
            }

            int count = incomingCounts[srcRank];
            if (count <= 0) {
                continue;
            }

            incomingByRank[srcRank].resize((size_t)count);

            MPI_Request req;
            MPI_Irecv(
                incomingByRank[srcRank].data(),
                count,
                mpiType,
                (int)srcRank,
                dataTag,
                MPI_COMM_WORLD,
                &req
            );
            payloadReceives.push_back(req);
        }

        vector<MPI_Request> payloadSends;

        for (unsigned int dstRank = 0; dstRank < mpiSize; dstRank++) {
            if (dstRank == mpiRank) {
                continue;
            }
            if (!willSendToRank[dstRank]) {
                continue;
            }

            int count = (int)toRank[dstRank].size();
            if (count <= 0) {
                continue;
            }

            MPI_Request req;
            MPI_Isend(
                toRank[dstRank].data(),
                count,
                mpiType,
                (int)dstRank,
                dataTag,
                MPI_COMM_WORLD,
                &req
            );
            payloadSends.push_back(req);
        }

        if (!payloadReceives.empty()) {
            MPI_Waitall((int)payloadReceives.size(), payloadReceives.data(), MPI_STATUSES_IGNORE);
        }
        if (!payloadSends.empty()) {
            MPI_Waitall((int)payloadSends.size(), payloadSends.data(), MPI_STATUSES_IGNORE);
        }


        // Flatten all received payloads
        received.clear();

        for (unsigned int srcRank = 0; srcRank < mpiSize; srcRank++) {
            vector<T>& list = incomingByRank[srcRank];
            for (unsigned int i = 0; i < list.size(); i++) {
                received.push_back(list[i]);
            }
        }

        // Clear outboxes
        for (unsigned int r = 0; r < mpiSize; r++) {
            toRank[r].clear();
        }
    }
};

#endif