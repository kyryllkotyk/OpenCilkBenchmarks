#include "pm2_bailinbailout.h"

/*<===========================================================================>

    Winter 2026
    Author: Kyryll Kotyk
    Advisor: Munehiro Fukuda

    This file includes various function that simplify phase implementation.
    The "glue" components

    Style Note:
        One parameter per line for function signatures
        Parantheses are closed on a new line, followed by body brace opening
        This rule does not apply if there are no parameters

<===========================================================================>*/

uint64_t BailInBailOut::percentFloorU64(
    uint64_t x,
    unsigned short pct
) {
    return (x * (uint64_t)pct) / 100ULL;
}

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
    unsigned int totalRanks
) {

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


unsigned int BailInBailOut::ownerRankFromGlobalId(
    unsigned int bankCountTotal,
    unsigned int totalRanks,
    unsigned int globalId
) {

    // Inverse mapping (global ID to owning rank, matches computeRange splitting)

    unsigned int base = bankCountTotal / totalRanks;
    unsigned int remainder = bankCountTotal % totalRanks;

    // Give +1 remainder to the first ranks until no remainder is left

    unsigned int cutoff = (base + 1) * remainder;

    if (globalId < cutoff) {
        return globalId / (base + 1);
    }
    else {
        return remainder + (globalId - cutoff) / base;
    }
}

