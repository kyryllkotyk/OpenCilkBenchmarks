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

unsigned int BailInBailOut::getInterestRate(
    uint64_t baseSeed,
    unsigned int run,
    unsigned int globalId,
    uint64_t minVal,
    uint64_t maxVal
) {
    uint64_t seed = makeSeed(
        baseSeed,
        run,
        0,
        globalId,
        STREAM_INTEREST_RATE_SELECTION
    );

    // Rebuild interest rates from seed to avoid MPI calls
    unsigned short rate = (unsigned short)randomInRangeFromSeed
    (seed, minVal, maxVal);

    return rate;
}

//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================
//=============================================================================

#include <sstream>
#include <iomanip>
#include <unordered_set>


void BailInBailOut::gatherAndPrintReportRank0Only(
    unsigned int mpiRank,
    unsigned int mpiSize,
    const std::string& localReport
) {
    int localLen = (int)localReport.size();
    vector<int> recvLens(mpiSize, 0);
    vector<int> displs(mpiSize, 0);

    MPI_Gather(
        &localLen,
        1,
        MPI_INT,
        recvLens.data(),
        1,
        MPI_INT,
        0,
        MPI_COMM_WORLD
    );

    int totalLen = 0;
    vector<char> allChars;

    if (mpiRank == 0) {
        for (unsigned int r = 0; r < mpiSize; r++) {
            displs[r] = totalLen;
            totalLen += recvLens[r];
        }
        allChars.resize((size_t)totalLen);
    }

    MPI_Gatherv(
        (void*)localReport.data(),
        localLen,
        MPI_CHAR,
        (mpiRank == 0 ? allChars.data() : nullptr),
        recvLens.data(),
        displs.data(),
        MPI_CHAR,
        0,
        MPI_COMM_WORLD
    );

    if (mpiRank == 0) {
        for (unsigned int r = 0; r < mpiSize; r++) {
            std::cout
                << "\n------------------------------ RANK " << r
                << " ------------------------------\n";
            std::cout.write(allChars.data() + displs[r], recvLens[r]);
            std::cout
                << "-----------------------------------------------------------------------\n";
        }
        std::cout.flush();
    }
}

uint64_t BailInBailOut::sumDebtU64(
    const vector<BailInBailOut::DebtEntry>& debts
) {
    uint64_t s = 0;
    for (const auto& e : debts) {
        s += e.amount;
    }
    return s;
}
void BailInBailOut::ASSERT_AND_REPORT_PRE_F(
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
) {
    std::ostringstream out;
    
    auto allreduceU64 = [](uint64_t localVal) -> uint64_t {
        uint64_t globalVal = 0;
        MPI_Allreduce(&localVal, &globalVal, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        return globalVal;
        };

    auto allreduceI64 = [](int64_t localVal) -> int64_t {
        int64_t globalVal = 0;
        MPI_Allreduce(&localVal, &globalVal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        return globalVal;
        };

    auto allreduceU64Max = [](uint64_t localVal) -> uint64_t {
        uint64_t globalVal = 0;
        MPI_Allreduce(&localVal, &globalVal, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        return globalVal;
        };

    auto allreduceI64Min = [](int64_t localVal) -> int64_t {
        int64_t globalVal = 0;
        MPI_Allreduce(&localVal, &globalVal, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
        return globalVal;
        };

    auto allreduceI64Max = [](int64_t localVal) -> int64_t {
        int64_t globalVal = 0;
        MPI_Allreduce(&localVal, &globalVal, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        return globalVal;
        };

    uint64_t localFailures = 0;
    uint64_t localWarnings = 0;

    auto fail = [&](const std::string& msg) {
        localFailures++;
        out << "FAIL: " << msg << "\n";
        };

    auto warn = [&](const std::string& msg) {
        localWarnings++;
        out << "WARN: " << msg << "\n";
        };

    out
        << "PRE-F STATE AUDIT (post-E3, pre-intervention)\n"
        << "run=" << run << " timestep=" << timestep
        << " | rank " << mpiRank << "/" << mpiSize << "\n"
        << "bankRange: start=" << bankGlobalStartIndex
        << " count=" << bankCountForRank
        << " | firmRange: start=" << firmGlobalStartIndex
        << " count=" << firmCountForRank
        << "\n\n";

    // -------------------------
    // 0) Shape checks (fatal)
    // -------------------------
    bool shapeOk = true;

    if (localBankLiquidity.size() != bankCountForRank) {
        fail("localBankLiquidity.size() != bankCountForRank");
        shapeOk = false;
    }
    if (localBankDebts.size() != bankCountForRank) {
        fail("localBankDebts.size() != bankCountForRank");
        shapeOk = false;
    }
    if (localFirmLiquidity.size() != firmCountForRank) {
        fail("localFirmLiquidity.size() != firmCountForRank");
        shapeOk = false;
    }
    if (localFirmDebts.size() != firmCountForRank) {
        fail("localFirmDebts.size() != firmCountForRank");
        shapeOk = false;
    }

    if (!localBankNeighbors.empty() && localBankNeighbors.size() != bankCountForRank) {
        warn("localBankNeighbors.size() != bankCountForRank (graph checks may be incomplete)");
    }

    // If shapes are broken, gather/print and optionally abort.
    if (!shapeOk) {
        uint64_t globalFailures = allreduceU64(localFailures);
        uint64_t globalWarnings = allreduceU64(localWarnings);

        if (mpiRank == 0) {
            std::cout
                << "================================================================================\n"
                << "PRE-F STATE AUDIT (post-E3, pre-intervention)\n"
                << "run=" << run << " timestep=" << timestep
                << " | ranks=" << mpiSize << "\n"
                << "RESULT: FAIL (shape mismatch)\n"
                << "failures=" << globalFailures << " warnings=" << globalWarnings << "\n"
                << "================================================================================\n";
        }

        gatherAndPrintReportRank0Only(mpiRank, mpiSize, out.str());

        if (hardAbortOnFailure && globalFailures != 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        return;
    }

    // -------------------------
    // 1) Scan bank + firm debts
    // -------------------------
    uint64_t localBankDebtEntries = 0;
    uint64_t localFirmDebtEntries = 0;

    uint64_t localZeroDebtEntries = 0;
    uint64_t localBadLenderIds = 0;
    uint64_t localBankSelfDebts = 0;

    uint64_t localNegativeLiquidityBanks = 0;
    uint64_t localInsolventBanks = 0; // liquidity < totalDebt

    int64_t localBankLiquiditySum = 0;
    int64_t localFirmLiquiditySum = 0;

    uint64_t localBankDebtSum = 0;
    uint64_t localFirmDebtSum = 0;

    int64_t localBankLiquidityMin = std::numeric_limits<int64_t>::max();
    int64_t localBankLiquidityMax = std::numeric_limits<int64_t>::min();


    struct WorstBank {
        unsigned int bankGlobalId;
        int64_t liquidity;
        uint64_t debtSum;
        int64_t deficit; // debtSum - liquidity
    };

    vector<WorstBank> localWorstByDeficit;
    localWorstByDeficit.reserve(32);

    // Banks
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;

        int64_t liq = localBankLiquidity[localB];
        localBankLiquiditySum += liq;

        if (liq < localBankLiquidityMin) localBankLiquidityMin = liq;
        if (liq > localBankLiquidityMax) localBankLiquidityMax = liq;

        if (liq < 0) localNegativeLiquidityBanks++;

        const auto& debts = localBankDebts[localB];

        uint64_t debtSum = 0;
        for (const auto& e : debts) {
            localBankDebtEntries++;
            debtSum += e.amount;

            if (e.amount == 0) localZeroDebtEntries++;
            if (e.lenderBankGlobalId >= bankCountTotal) localBadLenderIds++;
            if (e.lenderBankGlobalId == gB) localBankSelfDebts++;
        }

        localBankDebtSum += debtSum;

        if ((long double)liq < (long double)debtSum) {
            localInsolventBanks++;
        }

        // Track worst deficits (debtSum - liquidity)
        // Using signed deficit so negative means "fine".
        int64_t deficit = (int64_t)debtSum - liq;
        if (deficit > 0) {
            localWorstByDeficit.push_back({ gB, liq, debtSum, deficit });
        }
    }

    // Firms
    for (unsigned int localF = 0; localF < firmCountForRank; localF++) {
        localFirmLiquiditySum += localFirmLiquidity[localF];

        const auto& debts = localFirmDebts[localF];
        uint64_t debtSum = 0;

        for (const auto& e : debts) {
            localFirmDebtEntries++;
            debtSum += e.amount;

            if (e.amount == 0) localZeroDebtEntries++;
            if (e.lenderBankGlobalId >= bankCountTotal) localBadLenderIds++;
        }

        localFirmDebtSum += debtSum;
    }

    // Strong integrity checks:
    if (localZeroDebtEntries != 0) {
        fail("Found zero-amount debt entries (should be compacted away)");
    }
    if (localBadLenderIds != 0) {
        fail("Found lenderBankGlobalId out of range");
    }
    if (localBankSelfDebts != 0) {
        fail("Found bank self-debt entries (lender == borrower)");
    }

    // -------------------------
    // 2) Graph sanity checks
    // -------------------------
    uint64_t localNeighborBadIds = 0;
    uint64_t localNeighborSelfEdges = 0;
    uint64_t localNeighborDuplicates = 0;

    uint64_t localDegreeSum = 0;
    uint64_t localDegreeMin = std::numeric_limits<uint64_t>::max();
    uint64_t localDegreeMax = 0;

    if (!localBankNeighbors.empty()) {
        unsigned int neighborCountToCheck = (unsigned int)std::min(
            (size_t)bankCountForRank,
            localBankNeighbors.size()
        );

        for (unsigned int localB = 0; localB < neighborCountToCheck; localB++) {
            unsigned int gB = bankGlobalStartIndex + localB;
            const auto& nbrs = localBankNeighbors[localB];

            uint64_t deg = (uint64_t)nbrs.size();
            localDegreeSum += deg;
            if (deg < localDegreeMin) localDegreeMin = deg;
            if (deg > localDegreeMax) localDegreeMax = deg;

            std::unordered_set<unsigned int> seen;
            seen.reserve(nbrs.size() * 2 + 1);

            for (unsigned int n : nbrs) {
                if (n >= bankCountTotal) localNeighborBadIds++;
                if (n == gB) localNeighborSelfEdges++;
                if (!seen.insert(n).second) localNeighborDuplicates++;
            }
        }

        if (localNeighborBadIds != 0) fail("Graph neighbor list contains out-of-range IDs");
        if (localNeighborSelfEdges != 0) warn("Graph neighbor list contains self edges (may be ok if intended)");
        if (localNeighborDuplicates != 0) warn("Graph neighbor list contains duplicates (wastes sampling slots)");
    }

    // -------------------------
    // RNG values this timestep
    // -------------------------
    out << "\n[LOCAL RNG VALUES]\n";

    out << "Bank interest rates (bankGlobalId | interestRate | employeeWageCost):\n";
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;
        unsigned short rate = (unsigned short)getInterestRate(
            baseSeed, run, gB,
            minInterestRate, maxInterestRate
        );
        int64_t employeeWageCost = (int64_t)bankEmployeeWage
            * (int64_t)bankWorkerCount
            * (int64_t)wageConsumptionPercent / 100;
        out << "  bankGID=" << gB
            << " interestRate=" << rate
            << " employeeWageCost=" << employeeWageCost << "\n";
    }

    out << "Firm shocks, profit multipliers, workforce cost"
        << " (firmGlobalId | shock | profitMultiplier | workforceCost):\n";

    for (unsigned int localF = 0; localF < firmCountForRank; localF++) {
        unsigned int gF = firmGlobalStartIndex + localF;

        uint64_t shockSeed = makeSeed(
            baseSeed, run, timestep, gF, STREAM_SHOCK_GENERATION
        );
        short shock = (short)randomInRangeFromSeed(
            shockSeed, shockMultiplierMin, shockMultiplierMax
        );
        bool isNegative = (shockSeed & 1ULL);
        if (isNegative) shock *= -1;

        uint64_t profitSeed = makeSeed(
            baseSeed, run, timestep, gF, STREAM_PROFIT_MULTIPLIER_GENERATION
        );
        short profitMultiplier = (short)randomInRangeFromSeed(
            profitSeed, profitMultiplierMin, profitMultiplierMax
        );

        // Count workers assigned to this firm
        uint64_t workersForFirm = 0;
        for (unsigned int w = 0; w < workerCountTotal; w++) {
            unsigned int assignedFirm = selectForWorker(
                baseSeed, w, firmCountTotal, STREAM_WORKER_FIRM_ASSIGN
            );
            if (assignedFirm == gF) workersForFirm++;
        }
        uint64_t workforceCost = workersForFirm * (uint64_t)wage;

        out << "  firmGID=" << gF
            << " shock=" << shock
            << " profitMultiplier=" << profitMultiplier
            << " workers=" << workersForFirm
            << " workforceCost=" << workforceCost << "\n";
    }

    out << "Worker bank assignments (workerGlobalId | bankGlobalId | deposit):\n";
    uint64_t depositPerW = ((uint64_t)wage * (uint64_t)(100 - wageConsumptionPercent)) / 100ULL;
    for (unsigned int w = 0; w < workerCountTotal; w++) {
        unsigned int assignedBank = selectForWorker(
            baseSeed, w, bankCountTotal, STREAM_WORKER_BANK_ASSIGN
        );
        unsigned int assignedFirm = selectForWorker(
            baseSeed, w, firmCountTotal, STREAM_WORKER_FIRM_ASSIGN
        );
        out << "  workerGID=" << w
            << " bankGID=" << assignedBank
            << " firmGID=" << assignedFirm
            << " deposit=" << depositPerW << "\n";
    }


    // -------------------------
    // 3) Global reductions
    // -------------------------
    uint64_t globalFailures = allreduceU64(localFailures);
    uint64_t globalWarnings = allreduceU64(localWarnings);

    uint64_t globalBankDebtEntries = allreduceU64(localBankDebtEntries);
    uint64_t globalFirmDebtEntries = allreduceU64(localFirmDebtEntries);

    uint64_t globalBankDebtSum = allreduceU64(localBankDebtSum);
    uint64_t globalFirmDebtSum = allreduceU64(localFirmDebtSum);

    int64_t globalBankLiquiditySum = allreduceI64(localBankLiquiditySum);
    int64_t globalFirmLiquiditySum = allreduceI64(localFirmLiquiditySum);

    uint64_t globalNegativeLiquidityBanks = allreduceU64(localNegativeLiquidityBanks);
    uint64_t globalInsolventBanks = allreduceU64(localInsolventBanks);

    int64_t globalBankLiquidityMin = allreduceI64Min(localBankLiquidityMin);
    int64_t globalBankLiquidityMax = allreduceI64Max(localBankLiquidityMax);

    uint64_t globalNeighborBadIds = allreduceU64(localNeighborBadIds);
    uint64_t globalNeighborSelfEdges = allreduceU64(localNeighborSelfEdges);
    uint64_t globalNeighborDuplicates = allreduceU64(localNeighborDuplicates);

    uint64_t globalDegreeSum = allreduceU64(localDegreeSum);
    uint64_t globalDegreeMin = 0;
    uint64_t globalDegreeMax = 0;

    /*if (!localBankNeighbors.empty()) {
        uint64_t globalDegMinTmp = 0;
        uint64_t globalDegMaxTmp = 0;
        MPI_Allreduce(&localDegreeMin, &globalDegMinTmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&localDegreeMax, &globalDegMaxTmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        globalDegreeMin = globalDegMinTmp;
        globalDegreeMax = globalDegMaxTmp;
    }*/

    uint64_t localMin = (bankCountForRank > 0) ? localDegreeMin : std::numeric_limits<uint64_t>::max();
    uint64_t localMax = (bankCountForRank > 0) ? localDegreeMax : 0;
    uint64_t globalDegMinTmp = 0;
    uint64_t globalDegMaxTmp = 0;
    MPI_Allreduce(&localMin, &globalDegMinTmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&localMax, &globalDegMaxTmp, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    globalDegreeMin = globalDegMinTmp;
    globalDegreeMax = globalDegMaxTmp;


    // -------------------------
    // 4) Rank-local “top worst” list
    // -------------------------
    std::sort(
        localWorstByDeficit.begin(),
        localWorstByDeficit.end(),
        [](const WorstBank& a, const WorstBank& b) {
            if (a.deficit != b.deficit) return a.deficit > b.deficit;
            return a.bankGlobalId < b.bankGlobalId;
        }
    );

    unsigned int showWorst = (unsigned int)std::min<size_t>(10, localWorstByDeficit.size());

    out << "\n[LOCAL TOTALS]\n";
    out << "bankLiquiditySum=" << localBankLiquiditySum
        << " | firmLiquiditySum=" << localFirmLiquiditySum << "\n";
    out << "bankDebtEntries=" << localBankDebtEntries
        << " | firmDebtEntries=" << localFirmDebtEntries << "\n";
    out << "bankDebtSum=" << localBankDebtSum
        << " | firmDebtSum=" << localFirmDebtSum << "\n";
    out << "negativeLiquidityBanks=" << localNegativeLiquidityBanks
        << " | insolventBanks(liq < totalDebt)=" << localInsolventBanks << "\n";
    out << "bankLiquidityMin=" << localBankLiquidityMin
        << " | bankLiquidityMax=" << localBankLiquidityMax << "\n";

    if (!localBankNeighbors.empty()) {
        out << "\n[LOCAL GRAPH]\n";
        out << "degreeMin=" << localDegreeMin
            << " | degreeMax=" << localDegreeMax
            << " | degreeAvg~=" << (bankCountForRank ? (double)localDegreeSum / (double)bankCountForRank : 0.0) << "\n";
        out << "badNeighborIds=" << localNeighborBadIds
            << " | selfEdges=" << localNeighborSelfEdges
            << " | duplicates=" << localNeighborDuplicates << "\n";
    }


    if (showWorst > 0) {
        out << "\n[LOCAL WORST BANK DEFICITS] (top " << showWorst << ")\n";
        out << "bankGlobalId | liquidity | debtSum | deficit(debt-liq)\n";
        for (unsigned int i = 0; i < showWorst; i++) {
            const auto& w = localWorstByDeficit[i];
            out << std::setw(11) << w.bankGlobalId << " | "
                << std::setw(9) << w.liquidity << " | "
                << std::setw(7) << w.debtSum << " | "
                << std::setw(16) << w.deficit << "\n";
        }
    }

    if (printAllBanks) {
        out << "\n[LOCAL BANK LINES]\n";
        out << "bankGlobalId | liquidity | debtCount | debtSum | insolvent\n";
        for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
            unsigned int gB = bankGlobalStartIndex + localB;
            int64_t liq = localBankLiquidity[localB];
            uint64_t debtSum = sumDebtU64(localBankDebts[localB]);
            bool insolvent = ((long double)liq < (long double)debtSum);

            out << std::setw(11) << gB << " | "
                << std::setw(9) << liq << " | "
                << std::setw(8) << localBankDebts[localB].size() << " | "
                << std::setw(7) << debtSum << " | "
                << (insolvent ? "YES" : " no") << "\n";
        }
    }


    // -------------------------
    // 5) Global header (rank 0)
    // -------------------------
    if (mpiRank == 0) {
        std::cout
            << "================================================================================\n"
            << "PRE-F STATE AUDIT (post-E3, pre-intervention)\n"
            << "run=" << run << " timestep=" << timestep << " | ranks=" << mpiSize << "\n"
            << "banks: total=" << bankCountTotal << " | firms: total=" << firmCountTotal << "\n"
            << "RESULT: " << (globalFailures ? "FAIL" : "PASS")
            << " | failures=" << globalFailures
            << " warnings=" << globalWarnings << "\n"
            << "--------------------------------------------------------------------------------\n"
            << "[GLOBAL TOTALS]\n"
            << "bankLiquiditySum=" << globalBankLiquiditySum
            << " | firmLiquiditySum=" << globalFirmLiquiditySum << "\n"
            << "bankDebtEntries=" << globalBankDebtEntries
            << " | firmDebtEntries=" << globalFirmDebtEntries << "\n"
            << "bankDebtSum=" << globalBankDebtSum
            << " | firmDebtSum=" << globalFirmDebtSum << "\n"
            << "negativeLiquidityBanks=" << globalNegativeLiquidityBanks
            << " | insolventBanks(liq < totalDebt)=" << globalInsolventBanks << "\n"
            << "bankLiquidityMin=" << globalBankLiquidityMin
            << " | bankLiquidityMax=" << globalBankLiquidityMax << "\n";

        if (!localBankNeighbors.empty()) {
            std::cout
                << "--------------------------------------------------------------------------------\n"
                << "[GLOBAL GRAPH]\n"
                << "degreeMin=" << globalDegreeMin
                << " | degreeMax=" << globalDegreeMax
                << " | degreeAvg~=" << (bankCountTotal ? (double)globalDegreeSum / (double)bankCountTotal : 0.0) << "\n"
                << "badNeighborIds=" << globalNeighborBadIds
                << " | selfEdges=" << globalNeighborSelfEdges
                << " | duplicates=" << globalNeighborDuplicates << "\n";
        }

        std::cout
            << "================================================================================\n";
        std::cout.flush();
    }


    // Print rank blocks cleanly
    gatherAndPrintReportRank0Only(
        mpiRank,
        mpiSize,
        out.str()
    );

    if (hardAbortOnFailure && globalFailures != 0) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }


}