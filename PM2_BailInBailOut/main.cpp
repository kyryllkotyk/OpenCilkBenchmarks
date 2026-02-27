// main.cpp
#include <mpi.h>

#include <cstdio>
#include <cstdint>
#include <chrono>

// Change this include to your actual header
#include "pm2_bailinbailout.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int mpiRank = 0;
    int mpiSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    // ----------------------------
    // Small deterministic test case
    // ----------------------------
    // Keep these SMALL so you can print/debug if needed.
    // Use >1 rank to validate routing: try -n 2, -n 4.

    unsigned int runs = 1;
    unsigned int timesteps = 3;
    unsigned int logEvery = 1;
    uint64_t baseSeed = 123456789ULL;

    unsigned int bankCountTotal = 8;
    unsigned int firmCountTotal = 8;
    unsigned int workerCountTotal = 16;

    unsigned int initialBankLiquidity = 2000;
    unsigned int initialFirmLiquidity = 200;
    unsigned int initialProductionCost = 50;
    unsigned int wage = 1000;

    unsigned short wageConsumptionPercent = 30;
    unsigned short profitMultiplierMin = 95;
    unsigned short profitMultiplierMax = 105;
    unsigned short shockMultiplierMin = 0;
    unsigned short shockMultiplierMax = 5;

    unsigned short minInterestRate = 2;
    unsigned short maxInterestRate = 10;

    unsigned short interbankDensity = 25;
    unsigned short maxInterbankLenderSamplingK = 4;
    unsigned short maxInterbankLoanPercent = 30;

    unsigned short maxFirmLoanPercent = 30;
    unsigned short firmLenderDegree = 4;
    unsigned short firmRepayPercent = 25;
    unsigned short bankRepayPercent = 25;

    unsigned short bailInCoveragePercent = 50;
    unsigned short bailOutCoveragePercent = 50;

    if (mpiRank == 0) {
        printf("=== BailInBailOut smoke test ===\n");
        printf("mpiSize=%d banks=%u firms=%u workers=%u timesteps=%u runs=%u\n",
            mpiSize, bankCountTotal, firmCountTotal, workerCountTotal, timesteps, runs);
    }

    BailInBailOut sim;

    MPI_Barrier(MPI_COMM_WORLD);
    auto start = std::chrono::high_resolution_clock::now();

    sim.runBenchmarkInAndOut(
        runs,
        timesteps,
        logEvery,
        baseSeed,

        bankCountTotal,
        firmCountTotal,
        workerCountTotal,

        initialBankLiquidity,
        initialFirmLiquidity,
        initialProductionCost,
        wage,

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

        bailInCoveragePercent,
        bailOutCoveragePercent,
        true
    );

    MPI_Barrier(MPI_COMM_WORLD);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Elapsed time: " << duration.count() << " ms\n";

    if (mpiRank == 0) {
        printf("=== Done ===\n");
    }

    MPI_Finalize();


    return 0;
}