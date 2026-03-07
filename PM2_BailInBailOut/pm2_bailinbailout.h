#ifndef PM2_BAIL_IN_BAIL_OUT_
#define PM2_BAIL_IN_BAIL_OUT_ 

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
    /**
    * @brief Executes the Bail-In Bail-Out financial contagion benchmark
    *
    * @details This function runs the complete simulation of a financial system
    * containing banks, firms, and workers. The simulation is distributed across
    * MPI ranks, where each rank owns a part of the workers, firms, and banks.
    * The benchmark models firm production, worker wage deposits, 
    * bank lending, interbank loans, and financial distress policies
    *
    * For each run, two policy modes are executed:
    *
    * Bail-In (interbank debts are forgiven by the lending banks)
    * Bail-Out (external liquidity is injected into insolvent banks)
    *
    * Each run consists of multiple timesteps. During every timestep the following
    * phases are executed:
    *
    * A: Firm production, shocks, profit realization, and loan repayment
    * B: Worker consumption and wage deposits
    * C: Worker deposits and firm repayments transferred to banks
    * D: Firms request loans from banks and banks process requests
    * E: Banks repay interbank loans, accrue interest, and borrow if insolvent
    * F: Insolvency detection and application of bail-in or bail-out policy
    *
    * The benchmark is designed to stress distributed systems through:
    *
    * Irregular message passing (loan requests and acceptances)
    * Graph-based communication (interbank graph network)
    * Distributed state updates
    * Deterministic randomization 
    *
    * All randomness is deterministically seeded so that identical results are 
    * produced regardless of MPI system configuration.
    *
    * @param runs Number of independent simulation runs to perform
    * @param timesteps Number of timesteps within each run
    * @param logEvery Interval in timesteps at which debugging and verification
    * output may be produced when debug mode is enabled
    * @param baseSeed Base seed used for random number generation
    * The seed is mixed with run number, timestep, IDs, and stream tags for
    * reproducibility
    *
    * @param bankCountTotal Total number of banks in the simulation
    * @param firmCountTotal Total number of firms in the simulation
    * @param workerCountTotal Total number of workers in the simulation
    * @param bankWorkerCount Number of employees assigned to each bank.
    * Determines wage payments made by banks every timestep
    *
    * @param initialBankLiquidity Initial liquidity assigned to every bank
    * @param initialFirmLiquidity Initial liquidity assigned to every firm
    * @param initialProductionCost Initial production cost assigned to firms.
    * Production cost changes over time due to shocks and profit multipliers
    * @param wage Wage paid to each firm worker per timestep
    * @param bankEmployeeWage Wage paid to each bank employee per timestep
    *
    * @param wageConsumptionPercent Percentage of wages consumed by workers.
    * The remaining portion is deposited into worker bank accounts
    * @param profitMultiplierMin Minimum profit multiplier applied to firm
    * production output (in percent)
    * @param profitMultiplierMax Maximum profit multiplier applied to firm
    * production output (in percent)
    * @param shockMultiplierMin Minimum shock multiplier applied to firm
    * production cost (in percent)
    * @param shockMultiplierMax Maximum shock multiplier applied to firm
    * production cost (in percent)
    * @param minInterestRate Minimum interest rate applied to loans (in percent)
    * @param maxInterestRate Maximum interest rate applied to loans (in percent)
    *
    * @param interbankDensity Density of the interbank lending graph (percentage
    * of possible edges used when constructing neighbor lists)
    * @param maxInterbankLenderSamplingK Maximum number of lenders sampled when a
    * bank tries borrowing from neighbors in the interbank network
    * @param maxInterbankLoanPercent Maximum percentage of a bank's liquidity that
    * can be loaned to another bank during interbank lending
    * @param maxFirmLoanPercent Maximum percentage of a bank's liquidity that may
    * be lent to firms in a single timestep
    * @param firmLenderDegree Number of banks each firm is 
    * connected to in the firm-bank graph
    * @param firmRepayPercent Percentage of outstanding firm debt repaid each
    * timestep when liquidity permits
    * @param bankRepayPercent Percentage of outstanding interbank debt repaid
    * each timestep when liquidity permits
    *
    * @param bailInCoveragePercent Percentage of a bank’s debts that
    * will be forgiven during a bail-in intervention
    *
    * @param bailOutCoveragePercent Percentage of a bank’s deficit covered by
    * external liquidity injection during a bail-out intervention
    *
    * @param debug Enables detailed correctness verification and diagnostic
    * logging during simulation execution.
    */
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

    /**************************************************************************
    * This header section defines all functions implemented in the file       *
    * pm2_bailinbailout_rng.cpp                                               *
    **************************************************************************/

    //===============================STRUCTURES================================

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

    //===============================FUNCTIONS=================================

    /*
    * @brief Computes a SplitMix64 hash of the input value
    *
    * @param x Input value to be mixed
    * @return 64 bit hashed value
    */
    static uint64_t splitmix64Hash(
        uint64_t x
    );

    /*
    * @brief Deterministically selects an index for a worker
    *
    * @details
    * Uses a seeded random generator to assign a worker to an entity 
    * (bank or firm). The selection depends only on the base
    * seed, the worker's global ID, and the provided stream tag. 
    * The assignment must remain identical regardless of system configuration
    *
    * @param baseSeed Base random seed for the simulation.
    * @param workerGlobalId Global ID of the worker being assigned
    * @param countTotal Total number of possible targets to choose from
    * @param streamTag RNG stream identifier
    * @return Selected index between 0 and countTotal - 1
    */
    unsigned int selectForWorker(
        uint64_t baseSeed,
        unsigned int workerGlobalId,
        unsigned int countTotal,
        uint64_t streamTag = 1
    );

    /*
    * @brief Generates a deterministic mixed seed from multiple parameters
    *
    * @details 
    * Combines the base seed with run number, timestep, entity ID, and a stream
    * identifier to produce a unique seed. The resulting seed is then passed
    * through a SplitMix64 
    *
    * This is to make sure that all random values in the simulation are fully
    * independent of the system configuration
    *
    * @param baseSeed Base random seed of the simulation
    * @param run Current run index
    * @param timestep Current timestep index
    * @param globalId Global ID of the entity associated with the seed
    * @param streamTag RNG stream identifier 
    * @return Mixed 64 bit seed value
    */
    static uint64_t makeSeed(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int timestep,
        unsigned int globalId,
        uint64_t streamTag
    );

    /*
    * Generates a percentage value from 0 to 100 based on the mixed seed
    * @note The seed has to be mixed prior to calling this function!
    */
    unsigned short percent0to100FromMixedSeed(
        uint64_t seed
    );

    /*
    * Generates a random value between [minVal, maxVal] based on the mixed seed
    * @note The seed has to be mixed prior to calling this function!
    */
    uint64_t randomInRangeFromSeed(
        uint64_t seed,
        uint64_t minVal,
        uint64_t maxVal
    );

    /**************************************************************************
    * This header section defines all functions implemented in the file       *
    * pm2_bailinbailout_utils.cpp                                             *
    **************************************************************************/

    //===============================STRUCTURES================================
    
    //NOTE:: This struct is also used for firms. 
    struct DebtEntry {
        // Global bank ID
        unsigned int lenderBankGlobalId;
        // How much is owed
        uint64_t amount;
    };

    //================================FUNCTIONS================================

    // Returns floor(x * pct / 100) using integer arithmetic   
    uint64_t percentFloorU64(
        uint64_t x,
        unsigned short pct
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

    // Returns the MPI rank that owns the entity with the given global ID
    unsigned int ownerRankFromGlobalId(
        unsigned int countTotal,
        unsigned int totalRanks,
        unsigned int globalId
    );

    // Generates an interest rate for a bank
    // Used to avoid communication to receive interest rates
    // The interest rate gets regenerated using same parameters as 
    // it was originally generated from
    unsigned int getInterestRate(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int globalId,
        uint64_t minVal,
        uint64_t maxVal
    );


    // Collects report strings from all MPI ranks and prints them on rank 0
    static void gatherAndPrintReportRank0Only(
        unsigned int mpiRank,
        unsigned int mpiSize,
        const std::string& localReport
    );

    // Returns the total debt amount by adding all entries together
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

    /**************************************************************************
    * This header section defines all functions implemented in the file       *
    * pm2_bailinbailout_firms.cpp                                             *
    **************************************************************************/

    //================================STRUCTURES===============================

    struct FirmLoanRequest {
        unsigned int firmGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountRequested;
        unsigned short lenderInterestRate;
    };

    // Message used to send aggregated workforce cost contributions to the rank
    // that owns the target firm
    struct a0FirmWorkforceDeltaMessage {
        unsigned int firmGlobalId;
        uint64_t delta;
    };

    //================================DATATYPES================================

    MPI_Datatype a0MakeFirmWorkforceDeltaType();

    //================================FUNCTIONS================================

    // Builds the fixed bank neighbor list for a firm 
    // The neighbor set is generated and written to outNeighbors
    void buildFirmNeighbors(
        uint64_t baseSeed,
        unsigned int run,
        unsigned int firmGlobalId,
        unsigned int bankCountTotal,
        unsigned int firmLenderDegree,
        uint64_t stream,
        vector<unsigned int>& outNeighbors
    );


    // Repays firm's outstanding bank debts proportionally based on configured
    // repayment percentage and available liquidity
    // Records incoming repayments per lender bank
    void repayFirmLoansProRata(
        int64_t& firmLiquidity,
        vector<DebtEntry>& debts,
        unsigned short firmRepayPercent,
        vector<int64_t>& bankIncomingFromFirmRepay
    );

    // Computes the total workforce cost owed by each local firm
    // Routes combined wage deltas to the ranks that own those firms
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

    // Creates firm loan requests for firms with negative liquidity 
    // Routes each request to the rank that owns lender bank
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

    /**************************************************************************
    * This header section defines all functions implemented in the file       *
    * pm2_bailinbailout_banks.cpp                                             *
    **************************************************************************/

    //===============================STRUCTURES================================

    // Message used to apply a signed liquidity change to a bank
    // Used when routing worker deposits and firm repayments to bank owners.
    struct c1BankDeltaMessage {
        unsigned int bankGlobalId;
        int64_t delta;
    };

    // Acceptance message sent from lender bank back to firm that requested loan
    // All IDs are global
    struct d2FirmLoanAcceptance {
        unsigned int firmGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountGranted;
    };

    // Message representing an interbank loan repayment sent to the lender bank
    // The amount is credited to the bank
    struct e1InterbankRepaymentMessage {
        unsigned int lenderBankGlobalId;
        uint64_t amount;
    };

    // Request sent by bank to potential lender bank for interbank liquidity
    // Includes the requested amount and the lender's interest rate
    struct e3InterbankLoanRequest {
        unsigned int borrowerBankGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountRequested;
        unsigned short lenderInterestRate;
    };

    // Acceptance message sent by lender bank after granting interbank loan
    // Borrower applies granted amount and records the new debt locally
    struct e3InterbankLoanAcceptance {
        unsigned int borrowerBankGlobalId;
        unsigned int lenderBankGlobalId;
        uint64_t amountGranted;
    };

    // Message used to share a bank's current liquidity with neighboring banks
    // Allows borrowers to estimate lender capacity before asking to borrow
    struct e3NeighborLiquidityMessage {
        unsigned int bankGlobalId;
        int64_t liquidity;
    };

    // Local helper record representing one possible interbank lender
    // Stores lender ID, interest rate, and current liquidity for ranking
    struct e3LenderOption {
        unsigned int lenderGlobalId;
        unsigned short interestRate;
        int64_t liquidity;
    };

    //================================DATATYPES================================

    MPI_Datatype d1MakeFirmLoanRequestType();
    MPI_Datatype d4MakeFirmLoanAcceptanceType();
    MPI_Datatype c1MakeBankDeltaMessageType();
    MPI_Datatype e1MakeInterbankRepaymentType();
    MPI_Datatype e3MakeInterbankLoanRequestType();
    MPI_Datatype e3MakeInterbankLoanAcceptanceType();
    MPI_Datatype e3MakeNeighborLiquidityMessageType();

    //================================FUNCTIONS================================

    // Transfers worker deposits and firm loan repayments to owning bank ranks
    // and applies the resulting liquidity changes to local banks
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

    // Sends firm loan requests to the ranks that own the target banks and
    // collects all incoming requests that are processed locally
    void d1SendFirmLoanRequests(
        unsigned int mpiSize,
        vector<vector<FirmLoanRequest>>& firmLoanRequestsToRank,
        vector<FirmLoanRequest>& receivedLoanRequests
    );

    // Processes incoming firm loan requests for banks owned by this rank,
    // grants loans as possible (liquidity limits)
    void d2ProcessFirmLoanRequests(
        unsigned int firmCountTotal,
        unsigned int bankCountForRank,
        unsigned int bankGlobalStartIndex,
        unsigned int mpiSize,
        unsigned short maxFirmLoanPercent,
        vector<int64_t>& localBankLiquidity,
        vector<FirmLoanRequest>& receivedLoanRequests,
        vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank,
        vector<uint64_t>& bankFirmLoanOutflow
    );

    // Sends loan acceptance messages from lender banks back to the ranks
    // that own the borrowing firms
    void d4SendFirmLoanAcceptances(
        unsigned int mpiSize,
        vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank,
        vector<d2FirmLoanAcceptance>& receivedAcceptances
    );

    // Applies accepted firm loans by increasing firm liquidity and recording
    // the new debt entries 
    void d5ApplyFirmLoanAcceptances(
        unsigned int firmGlobalStartIndex,
        vector<int64_t>& localFirmLiquidity,
        vector<vector<DebtEntry>>& localFirmDebts,
        vector<d2FirmLoanAcceptance>& receivedAcceptances
    );

    // Repays outstandi ng interbank debts proportionally using available bank
    // liquidity and routes repayment messages to lender banks
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

    // Applies interest to all existing debts using the regenerated interest
    // rates per bank
    void e2ApplyInterestOnAllLoans(
        vector<vector<DebtEntry>>& eLocalBankDebts,
        vector<vector<DebtEntry>>& localFirmDebts,
        uint64_t baseSeed,
        unsigned int run,
        uint64_t minVal,
        uint64_t maxVal
    );

    // Banks with negative liquidity borrow from neighboring banks
    // in the interbank graph, up to sampling limits and other constraints
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

    /**************************************************************************
    * This header section defines and implements                              *
    * all templated functions used by the program                             *
    **************************************************************************/

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
