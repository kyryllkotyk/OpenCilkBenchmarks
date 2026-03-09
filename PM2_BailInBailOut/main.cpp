#include <cstdio>
#include <cstdint>
#include <chrono>
#include <iostream>

#include <mpi.h>
#include "pm2_bailinbailout.h"

static void printTimings(int mpiRank, int mpiSize, long long localMs) {
    std::vector<long long> allMs(mpiSize, 0);
    MPI_Gather(&localMs, 1, MPI_LONG_LONG,
        allMs.data(), 1, MPI_LONG_LONG,
        0, MPI_COMM_WORLD);

    if (mpiRank == 0) {
        std::vector<std::pair<long long, int>> ranked;
        ranked.reserve(mpiSize);
        for (int r = 0; r < mpiSize; r++)
            ranked.push_back({ allMs[r], r });
        std::sort(ranked.begin(), ranked.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

        printf("Timings (slowest first):\n");
        for (const auto& p : ranked)
            printf("  rank %d: %lld ms\n", p.second, p.first);
    }
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int mpiRank = 0;
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    unsigned int runs = 1;
    unsigned int timesteps = 8;
    unsigned int logEvery = 1;         
    uint64_t baseSeed = 123456789ULL;

    unsigned int bankCountTotal = 4;
    unsigned int firmCountTotal = 8;
    unsigned int workerCountTotal = 16;
    unsigned int bankWorkerCount = 2;

    unsigned int initialBankLiquidity = 5000;
    unsigned int initialFirmLiquidity = 100;
    unsigned int initialProductionCost = 200;
    unsigned int wage = 500;
    unsigned int bankEmployeeWage = 100;

    unsigned short wageConsumptionPercent = 50;   
    unsigned short profitMultiplierMin = 100;  
    unsigned short profitMultiplierMax = 100;
    unsigned short shockMultiplierMin = 0;   
    unsigned short shockMultiplierMax = 0;
    unsigned short minInterestRate = 7;    
    unsigned short maxInterestRate = 7;

    unsigned short interbankDensity = 0;   
    unsigned short maxInterbankLenderSamplingK = 1;
    unsigned short maxInterbankLoanPercent = 0;   
    unsigned short maxFirmLoanPercent = 20;
    unsigned short firmLenderDegree = 1;   
    unsigned short firmRepayPercent = 25;  
    unsigned short bankRepayPercent = 0;   
    unsigned int   interventionDelay = 1;
    unsigned short bailInCoveragePercent = 0;  
    unsigned short bailOutCoveragePercent = 0; 

    if (mpiRank == 0) {
        printf("=== BailInBailOut deterministic correctness test ===\n");
        printf("mpiSize=%d  banks=%u  firms=%u  workers=%u  "
            "timesteps=%u  logEvery=%u  runs=%u\n",
            mpiSize, bankCountTotal, firmCountTotal,
            workerCountTotal, timesteps, logEvery, runs);
        printf("\n--- Pre-computed expected globals (verify against JSONL) ---\n");
        uint64_t depPerW = (uint64_t)wage * (100 - wageConsumptionPercent) / 100ULL;
        uint64_t bankWageCost = (uint64_t)bankEmployeeWage * bankWorkerCount
            * wageConsumptionPercent / 100;
        printf("  depositPerWorker              = %llu\n", (unsigned long long)depPerW);
        printf("  bankEmployeeCost per bank/step= %llu\n", (unsigned long long)bankWageCost);
        printf("  total bank wage cost/step     = %llu\n", (unsigned long long)(bankWageCost * bankCountTotal));
        printf("  total worker deposits/step    = %llu\n", (unsigned long long)(depPerW * workerCountTotal));
        printf("  total firm revenue/step       = %llu\n", (unsigned long long)((uint64_t)firmCountTotal * initialProductionCost));
        printf("  total firm workforce cost     = %llu\n", (unsigned long long)((uint64_t)workerCountTotal * wage));
        printf("  initial bank liq sum          = %llu\n", (unsigned long long)((uint64_t)bankCountTotal * initialBankLiquidity));
        printf("  step-0 post_A expected sum    = %llu\n",
            (unsigned long long)((uint64_t)bankCountTotal * initialBankLiquidity
                - bankWageCost * bankCountTotal));
        printf("  step-0 post_C expected sum    = %llu (assuming repayments=0)\n",
            (unsigned long long)((uint64_t)bankCountTotal * initialBankLiquidity
                - bankWageCost * bankCountTotal
                + depPerW * workerCountTotal));
        fflush(stdout);
    }

    BailInBailOut sim;
    if (mpiRank == 0) {
        printf("=== Scenario 1: Baseline correctness (no interbank) ===\n");
        printf("  depositPerWorker              = 250\n");
        printf("  bankEmployeeCost per bank     = 100\n");
        printf("  step-0 expected snap_post_a   = 19600\n");
        printf("  step-0 expected snap_post_c   = 23600 (before firm repayments)\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto s1t0 = std::chrono::high_resolution_clock::now();

    sim.runBenchmarkInAndOut(
        runs,
        timesteps,
        logEvery,
        baseSeed,

        bankCountTotal,
        firmCountTotal,
        workerCountTotal,
        bankWorkerCount,

        initialBankLiquidity,
        initialFirmLiquidity,
        initialProductionCost,
        wage,
        bankEmployeeWage,

        wageConsumptionPercent,
        profitMultiplierMin,
        profitMultiplierMax,
        shockMultiplierMin,
        shockMultiplierMax,
        minInterestRate,
        maxInterestRate,

        interbankDensity,
        maxInterbankLenderSamplingK,
        maxInterbankLoanPercent,
        maxFirmLoanPercent,
        firmLenderDegree,
        firmRepayPercent,
        bankRepayPercent,

        interventionDelay,      

        bailInCoveragePercent,
        bailOutCoveragePercent,

        true  
    );

    auto s1t1 = std::chrono::high_resolution_clock::now();
    long long s1ms = std::chrono::duration_cast<std::chrono::milliseconds>(s1t1 - s1t0).count();

    if (mpiRank == 0) { 
        printf("Scenario 1 done\n");
        fflush(stdout); 
    }
    printTimings(mpiRank, mpiSize, s1ms);


    if (mpiRank == 0) {
        printf("\n=== Scenario 2: Interbank stress (E.1/E.2/E.3 active) ===\n");
        printf("  bankEmployeeCost per bank     = 200\n");
        printf("  initialBankLiquidity          = 50\n");
        printf("  After A.5.5 (before deposits) = -150 per bank\n");
        printf("  Bank with 0 workers depositing: -150   -> E.3 fires\n");
        printf("  Bank with 1 worker  depositing:  100   -> E.3 skips\n");
        printf("  Net debt factor per step: (1.07)*(0.75) = 0.8025 -> debt decays\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto s2t0 = std::chrono::high_resolution_clock::now();

    sim.runBenchmarkInAndOut(
        1,              // runs
        8,              // timesteps
        1,              // logEvery
        123456789ULL,   // baseSeed

        4,              // bankCountTotal
        4,              // firmCountTotal
        2,              // workerCountTotal  -> only 2 workers; 2 banks guaranteed 0 deposits
        2,              // bankWorkerCount

        50,             // initialBankLiquidity  -> low enough that A.5.5 goes negative
        100,            // initialFirmLiquidity
        200,            // initialProductionCost
        500,            // wage
        200,            // bankEmployeeWage  -> 200*2*50/100=200 cost/bank -> kills liquidity

        50,             // wageConsumptionPercent
        100, 100,       // profitMultiplier (pinned)
        0, 0,         // shockMultiplier  (pinned)
        7, 7,         // interestRate     (pinned)

        100,            // interbankDensity  -> full mesh (degree=3 for 4 banks)
        3,              // maxInterbankLenderSamplingK  -> sample all neighbors
        25,             // maxInterbankLoanPercent  -> lend 25% of liquidity
        20,             // maxFirmLoanPercent
        1,              // firmLenderDegree
        25,             // firmRepayPercent
        25,             // bankRepayPercent  -> E.1 repays 25% of interbank debt/step

        1,              // interventionDelay
        0,              // bailInCoveragePercent
        0,              // bailOutCoveragePercent

        true            // debug
    );

    auto s2t1 = std::chrono::high_resolution_clock::now();
    long long s2ms = std::chrono::duration_cast<std::chrono::milliseconds>(s2t1 - s2t0).count();

    if (mpiRank == 0) { printf("Scenario 2 complete.\n"); fflush(stdout); }
    printTimings(mpiRank, mpiSize, s2ms);

    if (mpiRank == 0)
        printf("\n=== Done ===\n");


    MPI_Finalize();
    return 0;
}
