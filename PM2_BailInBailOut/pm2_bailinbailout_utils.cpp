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
    unsigned int interventionDelay,

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

    if (firmLenderDegree < interventionDelay) {
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
    unsigned int countTotal,
    unsigned int totalRanks,
    unsigned int globalId
) {

    // Inverse mapping (global ID to owning rank, matches computeRange splitting)

    unsigned int base = countTotal / totalRanks;
    unsigned int remainder = countTotal % totalRanks;

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

/*<===========================================================================>

    JSONL_AUDIT_TIMESTEP
    ====================
    Exhaustive per-timestep JSONL snapshot covering both the post-E3/pre-F
    checkpoint and the post-F checkpoint in a single file.

    The caller must snapshot pre-F bank state BEFORE calling
    f1f2UpdateBankDistressAndApplyIntervention, then call this function
    AFTER it completes, passing both the snapshots and the live post-F state.

    Filename: BIBO_HHMMSSDDMMYYYY_TIMESTEP<6-padded>.jsonl

    Record types (one JSON object per line, sorted by record_type then gid):
        "meta"           -- simulation parameters + stream tag constants
        "bank_pre_f"     -- per-bank state BEFORE phase F
        "bank_post_f"    -- per-bank state AFTER phase F (debts/liq mutated)
        "firm"           -- per-firm state (unchanged by F)
        "worker"         -- per-worker assignment + deposit (local slice)
        "phase_f_summary"-- bail money, grace money, policy, distress counts
        "global"         -- global invariants, integrity violations, pass/fail

    Output is identical regardless of MPI rank count:
        - Each rank contributes only its owned slice of banks/firms/workers
        - Rank 0 gathers all lines via MPI_Gatherv
        - Lines are sorted by (record_type order, gid) before writing

<===========================================================================>*/

#include <ctime>
#include <fstream>
#include <unordered_set>

// ---------------------------------------------------------------------------
// Internal: escape a string for JSON
// ---------------------------------------------------------------------------
static std::string auditJsonEscape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 4);
    for (char c : s) {
        switch (c) {
        case '"':  out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\n': out += "\\n";  break;
        case '\r': out += "\\r";  break;
        case '\t': out += "\\t";  break;
        default:   out += c;      break;
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// Internal: build the per-run-policy subdirectory and return its path.
// Called once per (run, policy) pair in run.cpp before the timestep loop.
// Rank 0 creates the directory; path is broadcast so all ranks agree.
//
// Directory: DebugFiles/BIBO_HHMMSSmmm_DDMMYYYY_RUN<6>_BAILIN|BAILOUT/
// ---------------------------------------------------------------------------
#include <sys/stat.h>
#include <chrono>

std::string BailInBailOut::auditBuildRunFolder(
    unsigned int mpiRank,
    unsigned int mpiSize,
    unsigned int run,
    unsigned int policy
) {
    std::string folder;
    if (mpiRank == 0) {
        auto now = std::chrono::system_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;
        time_t t_raw = std::chrono::system_clock::to_time_t(now);
        struct tm* t = localtime(&t_raw);

        char buf[160];
        snprintf(buf, sizeof(buf),
            "DebugFiles/BIBO_%02d%02d%02d%03lld_%02d%02d%04d_RUN%06u_%s",
            t->tm_hour, t->tm_min, t->tm_sec,
            (long long)ms.count(),
            t->tm_mday, t->tm_mon + 1, t->tm_year + 1900,
            run,
            (policy == 0 ? "BAILIN" : "BAILOUT")
        );
        folder = std::string(buf);

        mkdir("DebugFiles", 0755);
        if (mkdir(folder.c_str(), 0755) != 0) {
            fprintf(stderr,
                "[AUDIT] Warning: could not create folder %s\n",
                folder.c_str());
        }
    }

    int len = (int)folder.size();
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    folder.resize((size_t)len);
    MPI_Bcast((void*)folder.data(), len, MPI_CHAR, 0, MPI_COMM_WORLD);
    return folder;
}

// ---------------------------------------------------------------------------
// Internal: build the full file path for one timestep inside the run folder.
// ---------------------------------------------------------------------------
static std::string auditBuildFilename(
    const std::string& runFolder,
    unsigned int timestep,
    unsigned int policy
) {
    char buf[64];
    snprintf(buf, sizeof(buf),
        "/%06u_%s.jsonl",
        timestep,
        (policy == 0 ? "BAILIN" : "BAILOUT")
    );
    return runFolder + std::string(buf);
}

// ---------------------------------------------------------------------------
// Internal: gather all rank JSONL blocks to rank 0 via MPI_Gatherv
// Returns the full concatenated string on rank 0, empty on other ranks
// ---------------------------------------------------------------------------
static std::string auditGatherToRank0(
    unsigned int mpiRank,
    unsigned int mpiSize,
    const std::string& localBlock
) {
    int localLen = (int)localBlock.size();
    std::vector<int> recvLens(mpiSize, 0);
    std::vector<int> displs(mpiSize, 0);

    MPI_Gather(
        &localLen, 1, MPI_INT,
        recvLens.data(), 1, MPI_INT,
        0, MPI_COMM_WORLD
    );

    std::string gathered;
    if (mpiRank == 0) {
        int total = 0;
        for (unsigned int r = 0; r < mpiSize; r++) {
            displs[r] = total;
            total += recvLens[r];
        }
        gathered.resize((size_t)total);
    }

    MPI_Gatherv(
        (void*)localBlock.data(), localLen, MPI_CHAR,
        (mpiRank == 0 ? (void*)gathered.data() : nullptr),
        recvLens.data(), displs.data(), MPI_CHAR,
        0, MPI_COMM_WORLD
    );

    return gathered;
}

// ---------------------------------------------------------------------------
// Internal: split gathered block into lines, sort deterministically, write
//
// Sort order (record_type):
//   0=meta, 1=bank_pre_f, 2=bank_post_f, 3=firm, 4=worker,
//   5=phase_f_summary, 6=global, 7=other
// Within same type: sorted by "gid" ascending (numeric)
// ---------------------------------------------------------------------------
static void auditSortAndWrite(
    std::ofstream& file,
    const std::string& allData
) {
    std::vector<std::string> lines;
    {
        size_t start = 0;
        while (start < allData.size()) {
            size_t end = allData.find('\n', start);
            if (end == std::string::npos) end = allData.size();
            if (end > start)
                lines.push_back(allData.substr(start, end - start));
            start = end + 1;
        }
    }

    auto extractStr = [](const std::string& line,
        const char* key) -> std::string {
            size_t p = line.find(key);
            if (p == std::string::npos) return "";
            p += strlen(key);
            size_t q = line.find('"', p);
            if (q == std::string::npos) return "";
            return line.substr(p, q - p);
        };

    auto extractGid = [](const std::string& line) -> uint64_t {
        const char* key = "\"gid\":";
        size_t p = line.find(key);
        if (p == std::string::npos) return UINT64_MAX;
        p += strlen(key);
        while (p < line.size() && line[p] == ' ') p++;
        uint64_t v = 0;
        while (p < line.size() && line[p] >= '0' && line[p] <= '9') {
            v = v * 10 + (uint64_t)(line[p] - '0');
            p++;
        }
        return v;
        };

    auto typeOrder = [](const std::string& t) -> int {
        if (t == "meta")            return 0;
        if (t == "bank_snap_post_a") return 1;
        if (t == "bank_snap_post_c") return 2;
        if (t == "bank_snap_post_d") return 3;
        if (t == "bank_pre_f")      return 4;
        if (t == "bank_post_f")     return 5;
        if (t == "firm")            return 6;
        if (t == "worker")          return 7;
        if (t == "phase_f_summary") return 8;
        if (t == "global")          return 9;
        return 10;
        };

    std::sort(lines.begin(), lines.end(),
        [&](const std::string& a, const std::string& b) {
            std::string ta = extractStr(a, "\"record_type\":\"");
            std::string tb = extractStr(b, "\"record_type\":\"");
            int oa = typeOrder(ta), ob = typeOrder(tb);
            if (oa != ob) return oa < ob;
            return extractGid(a) < extractGid(b);
        }
    );

    for (const auto& line : lines)
        file << line << '\n';
}

// ---------------------------------------------------------------------------
// Internal: emit one bank record (shared by pre-F and post-F paths)
// ---------------------------------------------------------------------------
void BailInBailOut::auditEmitBankRecord(
    std::ostringstream& out,
    const char* recordType,
    unsigned int gB,
    unsigned int mpiRank,
    int64_t liq,
    const std::vector<BailInBailOut::DebtEntry>& debts,
    unsigned int distressCount,
    const std::vector<unsigned int>& neighbors,
    unsigned int bankCountTotal,
    unsigned short interestRate,
    int64_t employeeWageCost,
    uint64_t seedInterest,
    uint64_t seedBankGraph,
    uint64_t xo0, uint64_t xo1, uint64_t xo2, uint64_t xo3
) {
    uint64_t debtSum = 0;
    for (const auto& e : debts) debtSum += e.amount;
    bool insolvent = (liq < (int64_t)debtSum);

    out << "{"
        << "\"record_type\":\"" << recordType << "\","
        << "\"gid\":" << gB << ","
        << "\"owner_rank\":" << mpiRank << ","
        << "\"liquidity\":" << liq << ","
        << "\"debt_count\":" << debts.size() << ","
        << "\"debt_sum\":" << debtSum << ","
        << "\"insolvent\":" << (insolvent ? "true" : "false") << ","
        << "\"distress_count\":" << distressCount << ","
        << "\"interest_rate\":" << interestRate << ","
        << "\"employee_wage_cost\":" << employeeWageCost << ","
        << "\"seed_interest_rate\":" << seedInterest << ","
        << "\"seed_bank_graph\":" << seedBankGraph << ","
        << "\"xo256_interest_seed_out0\":" << xo0 << ","
        << "\"xo256_interest_seed_out1\":" << xo1 << ","
        << "\"xo256_interest_seed_out2\":" << xo2 << ","
        << "\"xo256_interest_seed_out3\":" << xo3 << ","
        << "\"debts\":[";
    for (size_t d = 0; d < debts.size(); d++) {
        if (d > 0) out << ",";
        out << "{\"lender_gid\":" << debts[d].lenderBankGlobalId
            << ",\"amount\":" << debts[d].amount << "}";
    }
    out << "],"
        << "\"neighbors\":[";
    for (size_t n = 0; n < neighbors.size(); n++) {
        if (n > 0) out << ",";
        out << neighbors[n];
    }
    out << "]}\n";
}


void BailInBailOut::JSONL_AUDIT_TIMESTEP(
    unsigned int run,
    unsigned int timestep,
    unsigned int policy,
    unsigned int mpiRank,
    unsigned int mpiSize,
    const std::string& runFolder,

    unsigned int bankCountTotal,
    unsigned int bankGlobalStartIndex,
    unsigned int bankCountForRank,
    unsigned int firmCountTotal,
    unsigned int firmGlobalStartIndex,
    unsigned int firmCountForRank,
    unsigned int workerCountForRank,
    unsigned int workerGlobalStartIndex,

    // Mid-step phase snapshots (bank liquidity only - debts unchanged at A/C/D)
    const vector<int64_t>& snapPostABankLiq,
    const vector<int64_t>& snapPostCBankLiq,
    const vector<int64_t>& snapPostDBankLiq,

    // Pre-F snapshots (copied before f1f2 runs)
    const vector<int64_t>& preFBankLiquidity,
    const vector<vector<DebtEntry>>& preFBankDebts,
    const vector<unsigned int>& preFBankDistressCount,

    // Post-F live state
    vector<int64_t>& localBankLiquidity,
    vector<vector<DebtEntry>>& localBankDebts,
    vector<unsigned int>& postFBankDistressCount,

    // Phase F outputs
    uint64_t bailMoneyThisStep,
    uint64_t graceMoneyThisStep,

    // Firm state (unchanged by F)
    vector<int64_t>& localFirmLiquidity,
    vector<vector<DebtEntry>>& localFirmDebts,
    vector<vector<unsigned int>>& localBankNeighbors,
    vector<vector<unsigned int>>& localFirmNeighbors,
    vector<uint64_t>& localFirmProductionCost,

    // Simulation parameters
    uint64_t baseSeed,
    unsigned int bankWorkerCount,
    unsigned int bankEmployeeWage,
    unsigned int wage,
    unsigned int workerCountTotal,
    const vector<unsigned short>& localWorkerBankID,
    const vector<unsigned short>& localWorkerFirmID,
    unsigned short minInterestRate,
    unsigned short maxInterestRate,
    unsigned short shockMultiplierMin,
    unsigned short shockMultiplierMax,
    unsigned short profitMultiplierMin,
    unsigned short profitMultiplierMax,
    unsigned short wageConsumptionPercent,
    const vector<uint64_t>& localFirmWorkforceCost
) {
    // -----------------------------------------------------------------------
    // Build filename from pre-computed run folder (broadcast not needed here
    // since runFolder was already broadcast by auditBuildRunFolder in run.cpp)
    // -----------------------------------------------------------------------
    std::string filename = auditBuildFilename(runFolder, timestep, policy);

    std::ostringstream local;

    uint64_t depositPerWorker =
        ((uint64_t)wage * (uint64_t)(100 - wageConsumptionPercent)) / 100ULL;

    // ===========================
    // META (rank 0 only)
    // ===========================
    if (mpiRank == 0) {
        local
            << "{\"record_type\":\"meta\","
            << "\"gid\":0,"
            << "\"run\":" << run << ","
            << "\"timestep\":" << timestep << ","
            << "\"policy\":" << policy << ","
            << "\"policy_name\":\"" << (policy == 0 ? "bail-in" : "bail-out") << "\","
            << "\"mpi_ranks\":" << mpiSize << ","
            << "\"base_seed\":" << baseSeed << ","
            << "\"bank_count_total\":" << bankCountTotal << ","
            << "\"firm_count_total\":" << firmCountTotal << ","
            << "\"worker_count_total\":" << workerCountTotal << ","
            << "\"wage\":" << wage << ","
            << "\"bank_employee_wage\":" << bankEmployeeWage << ","
            << "\"bank_worker_count\":" << bankWorkerCount << ","
            << "\"wage_consumption_pct\":" << wageConsumptionPercent << ","
            << "\"deposit_per_worker\":" << depositPerWorker << ","
            << "\"min_interest_rate\":" << minInterestRate << ","
            << "\"max_interest_rate\":" << maxInterestRate << ","
            << "\"shock_min\":" << shockMultiplierMin << ","
            << "\"shock_max\":" << shockMultiplierMax << ","
            << "\"profit_min\":" << profitMultiplierMin << ","
            << "\"profit_max\":" << profitMultiplierMax << ","
            << "\"stream_bank_graph\":" << STREAM_BANK_GRAPH_BUILD << ","
            << "\"stream_worker_bank\":" << STREAM_WORKER_BANK_ASSIGN << ","
            << "\"stream_worker_firm\":" << STREAM_WORKER_FIRM_ASSIGN << ","
            << "\"stream_interest\":" << STREAM_INTEREST_RATE_SELECTION << ","
            << "\"stream_shock\":" << STREAM_SHOCK_GENERATION << ","
            << "\"stream_profit\":" << STREAM_PROFIT_MULTIPLIER_GENERATION << ","
            << "\"stream_firm_bank_graph\":" << STREAM_FIRM_BANK_GRAPH_BUILD
            << "}\n";
    }

    // ===========================
    // Per-bank constants (same for pre and post records, compute once)
    // ===========================
    struct BankConstants {
        unsigned short interestRate;
        int64_t        employeeWageCost;
        uint64_t       seedInterest;
        uint64_t       seedBankGraph;
        uint64_t       xo0, xo1, xo2, xo3;
    };

    std::vector<BankConstants> bankConst(bankCountForRank);
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;

        bankConst[localB].interestRate = (unsigned short)getInterestRate(
            baseSeed, run, gB, minInterestRate, maxInterestRate
        );
        bankConst[localB].employeeWageCost =
            (int64_t)bankEmployeeWage *
            (int64_t)bankWorkerCount *
            (int64_t)wageConsumptionPercent / 100;
        bankConst[localB].seedInterest = makeSeed(
            baseSeed, run, 0, gB, STREAM_INTEREST_RATE_SELECTION
        );
        bankConst[localB].seedBankGraph = makeSeed(
            baseSeed, run, 0, gB, STREAM_BANK_GRAPH_BUILD
        );

        Xoshiro256 rng(bankConst[localB].seedInterest);
        bankConst[localB].xo0 = rng.next();
        bankConst[localB].xo1 = rng.next();
        bankConst[localB].xo2 = rng.next();
        bankConst[localB].xo3 = rng.next();
    }

    // ===========================
    // BANK_SNAP_POST_A records (after A.5.5: employee wages deducted)
    // Invariant: liquidity == prev_post_f_liquidity - bankEmployeeCost
    // ===========================
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;
        local << "{"
            << "\"record_type\":\"bank_snap_post_a\","
            << "\"gid\":" << gB << ","
            << "\"owner_rank\":" << mpiRank << ","
            << "\"liquidity\":" << snapPostABankLiq[localB]
            << "}\n";
    }

    // ===========================
    // BANK_SNAP_POST_C records (after c1SendMoneyToBanks: deposits + repayments)
    // Invariant: sum == post_a_sum + total_deposits + total_firm_repayments
    // ===========================
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;
        local << "{"
            << "\"record_type\":\"bank_snap_post_c\","
            << "\"gid\":" << gB << ","
            << "\"owner_rank\":" << mpiRank << ","
            << "\"liquidity\":" << snapPostCBankLiq[localB]
            << "}\n";
    }

    // ===========================
    // BANK_SNAP_POST_D records (after D.3: firm loan outflows)
    // Invariant: sum == post_c_sum - total_loans_granted
    // E and F are no-ops with interbankDensity=0 and coverage=0,
    // so post_d_sum should equal post_f_bank_liquidity_sum
    // ===========================
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;
        local << "{"
            << "\"record_type\":\"bank_snap_post_d\","
            << "\"gid\":" << gB << ","
            << "\"owner_rank\":" << mpiRank << ","
            << "\"liquidity\":" << snapPostDBankLiq[localB]
            << "}\n";
    }

    // ===========================
    // BANK_PRE_F records
    // ===========================
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;
        const auto& bc = bankConst[localB];
        const std::vector<unsigned int>& nbrs =
            (localB < localBankNeighbors.size())
            ? localBankNeighbors[localB]
            : std::vector<unsigned int>{};

        auditEmitBankRecord(
            local, "bank_pre_f", gB, mpiRank,
            preFBankLiquidity[localB],
            preFBankDebts[localB],
            preFBankDistressCount[localB],
            nbrs, bankCountTotal,
            bc.interestRate, bc.employeeWageCost,
            bc.seedInterest, bc.seedBankGraph,
            bc.xo0, bc.xo1, bc.xo2, bc.xo3
        );
    }

    // ===========================
    // BANK_POST_F records
    // ===========================
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;
        const auto& bc = bankConst[localB];
        const std::vector<unsigned int>& nbrs =
            (localB < localBankNeighbors.size())
            ? localBankNeighbors[localB]
            : std::vector<unsigned int>{};

        auditEmitBankRecord(
            local, "bank_post_f", gB, mpiRank,
            localBankLiquidity[localB],
            localBankDebts[localB],
            postFBankDistressCount[localB],
            nbrs, bankCountTotal,
            bc.interestRate, bc.employeeWageCost,
            bc.seedInterest, bc.seedBankGraph,
            bc.xo0, bc.xo1, bc.xo2, bc.xo3
        );
    }

    // ===========================
    // FIRM records
    // ===========================
    for (unsigned int localF = 0; localF < firmCountForRank; localF++) {
        unsigned int gF = firmGlobalStartIndex + localF;

        int64_t  firmLiq = localFirmLiquidity[localF];
        uint64_t firmDebtSum = sumDebtU64(localFirmDebts[localF]);
        uint64_t prodCost = localFirmProductionCost[localF];
        uint64_t wfCost = localFirmWorkforceCost[localF];

        uint64_t shockSeed = makeSeed(
            baseSeed, run, timestep, gF, STREAM_SHOCK_GENERATION
        );
        uint64_t shockRaw = randomInRangeFromSeed(
            shockSeed, shockMultiplierMin, shockMultiplierMax
        );
        bool     shockNeg = (shockSeed & 1ULL);
        int64_t  shockEffective = shockNeg ? -(int64_t)shockRaw : (int64_t)shockRaw;

        uint64_t profitSeed = makeSeed(
            baseSeed, run, timestep, gF, STREAM_PROFIT_MULTIPLIER_GENERATION
        );
        uint64_t profitMultiplier = randomInRangeFromSeed(
            profitSeed, profitMultiplierMin, profitMultiplierMax
        );

        int64_t productionRevenue = (int64_t)(
            (double)prodCost * ((double)profitMultiplier / 100.0)
            );
        int64_t netPL = productionRevenue - (int64_t)wfCost;

        uint64_t seedFirmGraph = makeSeed(
            baseSeed, run, 0, gF, STREAM_FIRM_BANK_GRAPH_BUILD
        );

        Xoshiro256 rngF(shockSeed);
        uint64_t fxo0 = rngF.next(), fxo1 = rngF.next(),
            fxo2 = rngF.next(), fxo3 = rngF.next();

        const auto& fDebts = localFirmDebts[localF];

        local << "{"
            << "\"record_type\":\"firm\","
            << "\"gid\":" << gF << ","
            << "\"owner_rank\":" << mpiRank << ","
            << "\"liquidity\":" << firmLiq << ","
            << "\"debt_count\":" << fDebts.size() << ","
            << "\"debt_sum\":" << firmDebtSum << ","
            << "\"production_cost\":" << prodCost << ","
            << "\"workforce_cost\":" << wfCost << ","
            << "\"production_revenue\":" << productionRevenue << ","
            << "\"net_pl\":" << netPL << ","
            << "\"shock_seed\":" << shockSeed << ","
            << "\"shock_raw\":" << shockRaw << ","
            << "\"shock_negative\":" << (shockNeg ? "true" : "false") << ","
            << "\"shock_effective\":" << shockEffective << ","
            << "\"profit_seed\":" << profitSeed << ","
            << "\"profit_multiplier\":" << profitMultiplier << ","
            << "\"seed_firm_bank_graph\":" << seedFirmGraph << ","
            << "\"xo256_shock_seed_out0\":" << fxo0 << ","
            << "\"xo256_shock_seed_out1\":" << fxo1 << ","
            << "\"xo256_shock_seed_out2\":" << fxo2 << ","
            << "\"xo256_shock_seed_out3\":" << fxo3 << ","
            << "\"debts\":[";
        for (size_t d = 0; d < fDebts.size(); d++) {
            if (d > 0) local << ",";
            local << "{\"lender_gid\":" << fDebts[d].lenderBankGlobalId
                << ",\"amount\":" << fDebts[d].amount << "}";
        }
        local << "],"
            << "\"bank_neighbors\":[";
        if (localF < localFirmNeighbors.size()) {
            const auto& fnbrs = localFirmNeighbors[localF];
            for (size_t n = 0; n < fnbrs.size(); n++) {
                if (n > 0) local << ",";
                local << fnbrs[n];
            }
        }
        local << "]}\n";
    }

    // ===========================
    // WORKER records (local slice only - all ranks combined = full coverage)
    // ===========================
    for (unsigned int localW = 0; localW < workerCountForRank; localW++) {
        unsigned int gW = workerGlobalStartIndex + localW;

        unsigned int assignedBank = (unsigned int)localWorkerBankID[localW];
        unsigned int assignedFirm = (unsigned int)localWorkerFirmID[localW];

        uint64_t seedBank = makeSeed(baseSeed, 0, 0, gW, STREAM_WORKER_BANK_ASSIGN);
        uint64_t seedFirm = makeSeed(baseSeed, 0, 0, gW, STREAM_WORKER_FIRM_ASSIGN);

        local << "{"
            << "\"record_type\":\"worker\","
            << "\"gid\":" << gW << ","
            << "\"owner_rank\":" << mpiRank << ","
            << "\"assigned_bank_gid\":" << assignedBank << ","
            << "\"assigned_firm_gid\":" << assignedFirm << ","
            << "\"deposit_amount\":" << depositPerWorker << ","
            << "\"seed_bank_assign\":" << seedBank << ","
            << "\"seed_firm_assign\":" << seedFirm
            << "}\n";
    }

    // ===========================
    // PHASE_F_SUMMARY records (one per local bank, rank 0 gathers global)
    // ===========================
    // Emit per-bank F outcome: what changed between pre-F and post-F
    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;

        int64_t  preFliq = preFBankLiquidity[localB];
        int64_t  postFliq = localBankLiquidity[localB];
        uint64_t preDebt = sumDebtU64(preFBankDebts[localB]);
        uint64_t postDebt = sumDebtU64(localBankDebts[localB]);

        int64_t  liqDelta = postFliq - preFliq;
        int64_t  debtDelta = (int64_t)postDebt - (int64_t)preDebt;

        bool wasInsolvent = (preFliq < (int64_t)preDebt);
        bool stillInsolvent = (postFliq < (int64_t)postDebt);

        unsigned int preDistress = preFBankDistressCount[localB];
        unsigned int postDistress = postFBankDistressCount[localB];

        // Did F intervene this step for this bank?
        // Intervention happens when distress count >= interventionDelay
        // and distress % interventionDelay == 0.
        // We detect it by checking whether liq or debt actually changed
        bool intervened = (liqDelta != 0 || debtDelta != 0);

        local << "{"
            << "\"record_type\":\"phase_f_summary\","
            << "\"gid\":" << gB << ","
            << "\"owner_rank\":" << mpiRank << ","
            << "\"policy\":" << policy << ","
            << "\"policy_name\":\"" << (policy == 0 ? "bail-in" : "bail-out") << "\","
            << "\"pre_f_liquidity\":" << preFliq << ","
            << "\"post_f_liquidity\":" << postFliq << ","
            << "\"liquidity_delta\":" << liqDelta << ","
            << "\"pre_f_debt_sum\":" << preDebt << ","
            << "\"post_f_debt_sum\":" << postDebt << ","
            << "\"debt_delta\":" << debtDelta << ","
            << "\"was_insolvent\":" << (wasInsolvent ? "true" : "false") << ","
            << "\"still_insolvent\":" << (stillInsolvent ? "true" : "false") << ","
            << "\"pre_f_distress_count\":" << preDistress << ","
            << "\"post_f_distress_count\":" << postDistress << ","
            << "\"intervened\":" << (intervened ? "true" : "false")
            << "}\n";
    }

    // ===========================
    // GLOBAL invariants (reductions across all ranks, written by rank 0)
    // ===========================
    auto redU64 = [&](uint64_t v) -> uint64_t {
        uint64_t g = 0;
        MPI_Allreduce(&v, &g, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        return g;
        };
    auto redI64Sum = [&](int64_t v) -> int64_t {
        int64_t g = 0;
        MPI_Allreduce(&v, &g, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        return g;
        };
    auto redI64Min = [&](int64_t v) -> int64_t {
        int64_t g = 0;
        MPI_Allreduce(&v, &g, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
        return g;
        };
    auto redI64Max = [&](int64_t v) -> int64_t {
        int64_t g = 0;
        MPI_Allreduce(&v, &g, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        return g;
        };

    // --- Phase snapshot sums ---
    int64_t lSnapPostASum = 0, lSnapPostCSum = 0, lSnapPostDSum = 0;

    // --- Post-F bank stats ---
    int64_t  lBankLiqSum = 0, lBankLiqMin = INT64_MAX, lBankLiqMax = INT64_MIN;
    uint64_t lBankDebtSum = 0, lBankDebtCnt = 0;
    uint64_t lNegBanks = 0, lInsolventBanks = 0;
    uint64_t lZeroDebts = 0, lBadLenders = 0, lSelfDebts = 0;
    uint64_t lGraphBad = 0, lGraphSelf = 0, lGraphDupes = 0, lDegSum = 0;

    // --- Pre-F bank stats (for diff) ---
    int64_t  lPreBankLiqSum = 0;
    uint64_t lPreBankDebtSum = 0;
    uint64_t lPreInsolvent = 0;

    // --- Phase F totals ---
    uint64_t lBailMoney = bailMoneyThisStep;
    uint64_t lGraceMoney = graceMoneyThisStep;

    // --- Firm stats ---
    int64_t  lFirmLiqSum = 0, lFirmLiqMin = INT64_MAX, lFirmLiqMax = INT64_MIN;
    uint64_t lFirmDebtSum = 0, lFirmDebtCnt = 0;

    for (unsigned int localB = 0; localB < bankCountForRank; localB++) {
        unsigned int gB = bankGlobalStartIndex + localB;

        // Phase snapshot sums
        lSnapPostASum += snapPostABankLiq[localB];
        lSnapPostCSum += snapPostCBankLiq[localB];
        lSnapPostDSum += snapPostDBankLiq[localB];

        // Post-F
        int64_t liq = localBankLiquidity[localB];
        lBankLiqSum += liq;
        if (liq < lBankLiqMin) lBankLiqMin = liq;
        if (liq > lBankLiqMax) lBankLiqMax = liq;
        if (liq < 0) lNegBanks++;

        uint64_t ds = 0;
        for (const auto& e : localBankDebts[localB]) {
            lBankDebtCnt++;
            ds += e.amount;
            if (e.amount == 0)                         lZeroDebts++;
            if (e.lenderBankGlobalId >= bankCountTotal) lBadLenders++;
            if (e.lenderBankGlobalId == gB)             lSelfDebts++;
        }
        lBankDebtSum += ds;
        if (liq < (int64_t)ds) lInsolventBanks++;

        // Pre-F
        lPreBankLiqSum += preFBankLiquidity[localB];
        lPreBankDebtSum += sumDebtU64(preFBankDebts[localB]);
        if (preFBankLiquidity[localB] < (int64_t)sumDebtU64(preFBankDebts[localB]))
            lPreInsolvent++;

        // Graph
        if (localB < localBankNeighbors.size()) {
            std::unordered_set<unsigned int> seen;
            for (unsigned int n : localBankNeighbors[localB]) {
                lDegSum++;
                if (n >= bankCountTotal) lGraphBad++;
                if (n == gB)             lGraphSelf++;
                if (!seen.insert(n).second) lGraphDupes++;
            }
        }
    }

    for (unsigned int localF = 0; localF < firmCountForRank; localF++) {
        int64_t liq = localFirmLiquidity[localF];
        lFirmLiqSum += liq;
        if (liq < lFirmLiqMin) lFirmLiqMin = liq;
        if (liq > lFirmLiqMax) lFirmLiqMax = liq;
        for (const auto& e : localFirmDebts[localF]) {
            lFirmDebtCnt++;
            lFirmDebtSum += e.amount;
            if (e.amount == 0)                         lZeroDebts++;
            if (e.lenderBankGlobalId >= bankCountTotal) lBadLenders++;
        }
    }

    if (bankCountForRank == 0) { lBankLiqMin = 0; lBankLiqMax = 0; }
    if (firmCountForRank == 0) { lFirmLiqMin = 0; lFirmLiqMax = 0; }

    // Reductions
    int64_t  gBankLiqSum = redI64Sum(lBankLiqSum);
    int64_t  gBankLiqMin = redI64Min(lBankLiqMin);
    int64_t  gBankLiqMax = redI64Max(lBankLiqMax);
    uint64_t gBankDebtSum = redU64(lBankDebtSum);
    uint64_t gBankDebtCnt = redU64(lBankDebtCnt);
    uint64_t gNegBanks = redU64(lNegBanks);
    uint64_t gInsolventBanks = redU64(lInsolventBanks);
    uint64_t gZeroDebts = redU64(lZeroDebts);
    uint64_t gBadLenders = redU64(lBadLenders);
    uint64_t gSelfDebts = redU64(lSelfDebts);
    uint64_t gGraphBad = redU64(lGraphBad);
    uint64_t gGraphSelf = redU64(lGraphSelf);
    uint64_t gGraphDupes = redU64(lGraphDupes);
    uint64_t gDegSum = redU64(lDegSum);
    int64_t  gPreBankLiqSum = redI64Sum(lPreBankLiqSum);
    uint64_t gPreBankDebtSum = redU64(lPreBankDebtSum);
    uint64_t gPreInsolvent = redU64(lPreInsolvent);
    uint64_t gBailMoney = redU64(lBailMoney);
    uint64_t gGraceMoney = redU64(lGraceMoney);
    int64_t  gFirmLiqSum = redI64Sum(lFirmLiqSum);
    int64_t  gFirmLiqMin = redI64Min(lFirmLiqMin);
    int64_t  gFirmLiqMax = redI64Max(lFirmLiqMax);
    uint64_t gFirmDebtSum = redU64(lFirmDebtSum);
    uint64_t gFirmDebtCnt = redU64(lFirmDebtCnt);
    // Phase snapshot sums (global totals for step-by-step arithmetic verification)
    int64_t  gSnapPostASum = redI64Sum(lSnapPostASum);
    int64_t  gSnapPostCSum = redI64Sum(lSnapPostCSum);
    int64_t  gSnapPostDSum = redI64Sum(lSnapPostDSum);

    bool pass = (gZeroDebts == 0) && (gBadLenders == 0) && (gSelfDebts == 0)
        && (gGraphBad == 0);

    if (mpiRank == 0) {
        double avgDeg = bankCountTotal
            ? (double)gDegSum / (double)bankCountTotal : 0.0;

        local << "{"
            << "\"record_type\":\"global\","
            << "\"gid\":0,"
            << "\"run\":" << run << ","
            << "\"timestep\":" << timestep << ","
            << "\"pass\":" << (pass ? "true" : "false") << ","
            // Phase snapshot sums (verify step-by-step arithmetic)
            // post_a  = prev_post_f - bank_employee_cost_total
            // post_c  = post_a + total_deposits + total_firm_repayments
            // post_d  = post_c - total_loans_granted
            // post_f  = post_d  (when E and F are no-ops)
            << "\"snap_post_a_bank_liq_sum\":" << gSnapPostASum << ","
            << "\"snap_post_c_bank_liq_sum\":" << gSnapPostCSum << ","
            << "\"snap_post_d_bank_liq_sum\":" << gSnapPostDSum << ","
            // Post-F bank liquidity
            << "\"post_f_bank_liquidity_sum\":" << gBankLiqSum << ","
            << "\"post_f_bank_liquidity_min\":" << gBankLiqMin << ","
            << "\"post_f_bank_liquidity_max\":" << gBankLiqMax << ","
            << "\"post_f_banks_negative_liquidity\":" << gNegBanks << ","
            << "\"post_f_banks_insolvent\":" << gInsolventBanks << ","
            // Post-F bank debt
            << "\"post_f_bank_debt_entry_count\":" << gBankDebtCnt << ","
            << "\"post_f_bank_debt_sum\":" << gBankDebtSum << ","
            // Pre-F bank totals (for comparison)
            << "\"pre_f_bank_liquidity_sum\":" << gPreBankLiqSum << ","
            << "\"pre_f_bank_debt_sum\":" << gPreBankDebtSum << ","
            << "\"pre_f_banks_insolvent\":" << gPreInsolvent << ","
            // Phase F money flows
            << "\"bail_money_this_step\":" << gBailMoney << ","
            << "\"grace_money_this_step\":" << gGraceMoney << ","
            // Firm state
            << "\"firm_liquidity_sum\":" << gFirmLiqSum << ","
            << "\"firm_liquidity_min\":" << gFirmLiqMin << ","
            << "\"firm_liquidity_max\":" << gFirmLiqMax << ","
            << "\"firm_debt_entry_count\":" << gFirmDebtCnt << ","
            << "\"firm_debt_sum\":" << gFirmDebtSum << ","
            // Integrity violations
            << "\"violations_zero_debt_entries\":" << gZeroDebts << ","
            << "\"violations_bad_lender_ids\":" << gBadLenders << ","
            << "\"violations_bank_self_debts\":" << gSelfDebts << ","
            << "\"violations_graph_bad_ids\":" << gGraphBad << ","
            << "\"violations_graph_self_edges\":" << gGraphSelf << ","
            << "\"violations_graph_duplicates\":" << gGraphDupes << ","
            // Graph
            << "\"interbank_graph_avg_degree\":" << avgDeg << ","
            << "\"interbank_graph_total_edges\":" << gDegSum
            << "}\n";
    }

    // -----------------------------------------------------------------------
    // Gather all lines to rank 0, sort by (type, gid), write file
    // -----------------------------------------------------------------------
    std::string gathered = auditGatherToRank0(mpiRank, mpiSize, local.str());

    if (mpiRank == 0) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            fprintf(stderr,
                "JSONL_AUDIT_TIMESTEP: could not open %s for writing\n",
                filename.c_str());
            return;
        }
        auditSortAndWrite(file, gathered);
        file.close();
        printf("[AUDIT] Written: %s\n", filename.c_str());
        fflush(stdout);
    }
}