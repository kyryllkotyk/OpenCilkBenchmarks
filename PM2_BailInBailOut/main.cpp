// main.cpp
#include <mpi.h>

#include <cstdio>
#include <cstdint>
#include <chrono>

#include "pm2_bailinbailout.h"

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

    unsigned short bailInCoveragePercent = 0;
    unsigned short bailOutCoveragePercent = 0;


    if (mpiRank == 0) {
        printf("=== BailInBailOut smoke test ===\n");
        fflush(stdout);
        printf("mpiSize=%d banks=%u firms=%u workers=%u timesteps=%u runs=%u\n",
            mpiSize, bankCountTotal, firmCountTotal, workerCountTotal, timesteps, runs);
        fflush(stdout);
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

        bailInCoveragePercent,
        bailOutCoveragePercent,
        true
    );

    MPI_Barrier(MPI_COMM_WORLD);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    if (mpiRank == 0) {
        std::cout << "Elapsed time: " << duration.count() << " ms\n";
        printf("=== Done ===\n");
    }

    MPI_Finalize();


    return 0;
}
