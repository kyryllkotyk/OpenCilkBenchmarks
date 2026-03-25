#include "pm2_bailinbailout.h"
#include <unordered_map>

/*<===========================================================================>

    Winter 2026
    Author: Kyryll Kotyk
    Advisor: Munehiro Fukuda

    This file contains all Firm related logic, including lender sampling,
    requesting loans, repaying them, and building neighbor graph

    Style Note:
        One parameter per line for function signatures
        Parantheses are closed on a new line, followed by body brace opening
        This rule does not apply if there are no parameters

<===========================================================================>*/

void BailInBailOut::buildFirmNeighbors(
    uint64_t baseSeed,
    unsigned int run,
    unsigned int firmGlobalId,
    unsigned int bankCountTotal,
    unsigned int firmLenderDegree,
    uint64_t stream,
    vector<unsigned int>& outNeighbors
) {
    outNeighbors.clear();

    if (bankCountTotal == 0 || firmLenderDegree == 0) {
        return;
    }

    // Degree can't go above bankCountTotal
    unsigned int degree = min(firmLenderDegree, bankCountTotal);

    // Get the seed
    uint64_t seedy = makeSeed(
        baseSeed,
        run,
        0,
        firmGlobalId,
        stream
    );

    outNeighbors.reserve(degree);

    while (outNeighbors.size() < degree) {
        seedy = splitmix64Hash(seedy);

        unsigned int candidate = (unsigned int)(seedy % bankCountTotal);

        // Enforce uniqueness
        bool already = false;
        for (unsigned int x : outNeighbors) {
            if (x == candidate) {
                already = true;
                break;
            }
        }
        if (!already) {
            outNeighbors.push_back(candidate);
        }
    }
}


// Pro-rata repayment for one firm.
// Writes lender receipts into bankIncomingFromFirmRepay for its global index
void BailInBailOut::repayFirmLoansProRata(
    int64_t& firmLiquidity,
    vector<DebtEntry>& debts,
    unsigned short firmRepayPercent,
    vector<int64_t>& bankIncomingFromFirmRepay
) {
    // Return if not enough liquidity, no debts, or simulation set to not repay
    if (firmLiquidity <= 0 || debts.empty() || firmRepayPercent == 0) {
        return;
    }

    // How many debts are owed
    const size_t n = debts.size();
    vector<uint64_t> scheduled(n, 0);

    uint64_t scheduledTotal = 0;

    // Go through all debts
    for (size_t i = 0; i < n; i++) {
        // How much this specific debt is
        uint64_t debtAmount = debts[i].amount;
        // If this debt is zero, return
        if (debtAmount == 0) {
            continue;
        }

        uint64_t s = percentFloorU64(debtAmount, firmRepayPercent);
        scheduled[i] = s;
        scheduledTotal += s;
    }

    // Nothing to pay
    if (scheduledTotal == 0) {
        return;
    }

    // Maximum money the firm can use to pay
    uint64_t canPay = (uint64_t)firmLiquidity;
    // How much they are allowed to pay total
    uint64_t payTotal = (scheduledTotal <= canPay) ? scheduledTotal : canPay;

    if (payTotal == 0) {
        return;
    }


    vector<uint64_t> paymentReceived(n, 0);
    vector<uint64_t> proportionalRemainder(n, 0);

    // Tracker of how much got assigned so far
    uint64_t sumPay = 0;
    for (size_t i = 0; i < n; i++) {
        if (scheduled[i] == 0) {
            continue;
        }

        // Payment is the floor of (payTotal * scheduled / scheduledTotal)
        unsigned __int128 scaledNumerator = (unsigned __int128)payTotal
            * (unsigned __int128)scheduled[i];
        uint64_t paymentAmount = (uint64_t)(scaledNumerator / scheduledTotal);
        uint64_t proportionalRemain = (uint64_t)(scaledNumerator %
            scheduledTotal);

        // Safety cap
        if (paymentAmount > debts[i].amount) {
            paymentAmount = debts[i].amount;
        }

        paymentReceived[i] = paymentAmount;
        proportionalRemainder[i] = proportionalRemain;
        sumPay += paymentAmount;
    }

    // Distribute leftover 1 by 1 by largest remainder 
    uint64_t leftover = payTotal - sumPay;
    if (leftover > 0) {
        vector<size_t> order;
        order.reserve(n);

        for (size_t i = 0; i < n; i++) {
            if (scheduled[i] == 0) {
                continue;
            }
            if (debts[i].amount > paymentReceived[i]) {
                order.push_back(i);
            }
        }

        sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            if (proportionalRemainder[a] != proportionalRemainder[b]) {
                return proportionalRemainder[a] > proportionalRemainder[b];
            }
            if (debts[a].lenderBankGlobalId != debts[b].lenderBankGlobalId) {
                return debts[a].lenderBankGlobalId < debts[b].lenderBankGlobalId;
            }
            return a < b;
            });

        for (size_t k = 0; k < order.size() && leftover > 0; k++) {
            size_t i = order[k];
            if (debts[i].amount > paymentReceived[i]) {
                paymentReceived[i] += 1;
                leftover -= 1;
            }
        }
    }

    // Apply payments
    uint64_t actuallyPaid = 0;
    for (size_t i = 0; i < n; i++) {
        uint64_t p = paymentReceived[i];
        if (p == 0) {
            continue;
        }

        debts[i].amount -= p;
        bankIncomingFromFirmRepay[debts[i].lenderBankGlobalId] += (int64_t)p;
        actuallyPaid += p;
    }

    firmLiquidity -= (int64_t)actuallyPaid;

    // Cleanup paid debts
    for (size_t i = 0; i < debts.size(); ) {
        if (debts[i].amount == 0) {
            debts[i] = debts.back();
            debts.pop_back();
        }
        else {
            i++;
        }
    }
}

void BailInBailOut::a0BuildFirmPhaseAParameters(
    const uint64_t baseSeed,
    const unsigned int run,
    const unsigned int step,
    const unsigned int firmGlobalStartIndex,
    const unsigned int firmCountForRank,

    const unsigned short firmShockMinPercent,
    const unsigned short firmShockMaxPercent,
    const unsigned short firmProfitMinPercent,
    const unsigned short firmProfitMaxPercent,
    const unsigned short firmOperatingCostMinPercent,
    const unsigned short firmOperatingCostMaxPercent,
    const unsigned short firmShockProfitMultiplierMinPercent,
    const unsigned short firmShockProfitMultiplierMaxPercent,
    const unsigned short firmRareEventImpactMinPercent,
    const unsigned short firmRareEventImpactMaxPercent,

    const uint64_t streamFirmShockPercent,
    const uint64_t streamFirmProfitPercent,
    const uint64_t streamFirmOperatingCostPercent,
    const uint64_t streamFirmShockProfitMultiplierPercent,
    const uint64_t streamFirmRareEventProbability,
    const uint64_t streamFirmRareEventImpactPercent,
    const uint64_t streamFirmRareEventSign,

    vector<unsigned short>& localFirmShockMinPercent,
    vector<unsigned short>& localFirmShockMaxPercent,
    vector<unsigned short>& localFirmProfitMinPercent,
    vector<unsigned short>& localFirmProfitMaxPercent,
    vector<unsigned short>& localFirmOperatingCostPercent,
    vector<unsigned short>& localFirmShockProfitMultiplierPercent,
    vector<unsigned short>& localFirmRareEventRollPercent,
    vector<unsigned short>& localFirmRareEventImpactPercent,
    vector<bool>& localFirmRareEventIsPositive
) {
    for (unsigned int localFirm = 0; localFirm < firmCountForRank; localFirm++) {
        unsigned int firmGlobalId = firmGlobalStartIndex + localFirm;

        unsigned short shockA = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmShockPercent,
            firmShockMinPercent,
            firmShockMaxPercent
        );

        unsigned short shockB = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmShockPercent + 1ULL,
            firmShockMinPercent,
            firmShockMaxPercent
        );

        if (shockA <= shockB) {
            localFirmShockMinPercent[localFirm] = shockA;
            localFirmShockMaxPercent[localFirm] = shockB;
        }
        else {
            localFirmShockMinPercent[localFirm] = shockB;
            localFirmShockMaxPercent[localFirm] = shockA;
        }

        unsigned short profitA = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmProfitPercent,
            firmProfitMinPercent,
            firmProfitMaxPercent
        );

        unsigned short profitB = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmProfitPercent + 1ULL,
            firmProfitMinPercent,
            firmProfitMaxPercent
        );

        if (profitA <= profitB) {
            localFirmProfitMinPercent[localFirm] = profitA;
            localFirmProfitMaxPercent[localFirm] = profitB;
        }
        else {
            localFirmProfitMinPercent[localFirm] = profitB;
            localFirmProfitMaxPercent[localFirm] = profitA;
        }

        localFirmOperatingCostPercent[localFirm] = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmOperatingCostPercent,
            firmOperatingCostMinPercent,
            firmOperatingCostMaxPercent
        );

        localFirmShockProfitMultiplierPercent[localFirm] = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmShockProfitMultiplierPercent,
            firmShockProfitMultiplierMinPercent,
            firmShockProfitMultiplierMaxPercent
        );

        localFirmInitialLiquidityMultiplierPercent[localFirm] = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmInitialLiquidityMultiplierPercent,
            firmInitialLiquidityMultiplierMinPercent,
            firmInitialLiquidityMultiplierMaxPercent
        );

        localFirmRareEventRollPercent[localFirm] = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmRareEventProbability,
            0,
            100
        );

        localFirmRareEventImpactPercent[localFirm] = randomPercentFromStream(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmRareEventImpactPercent,
            firmRareEventImpactMinPercent,
            firmRareEventImpactMaxPercent
        );

        uint64_t rareEventSignSeed = makeSeed(
            baseSeed,
            run,
            step,
            firmGlobalId,
            streamFirmRareEventSign
        );

        localFirmRareEventIsPositive[localFirm] = ((rareEventSignSeed & 1ULL) != 0ULL);

    }
}

MPI_Datatype BailInBailOut::a0MakeFirmWorkforceDeltaType(){
    MPI_Datatype t;

    int blockLengths[2] = { 1, 1 };
    MPI_Aint offsets[2];
    offsets[0] = (MPI_Aint)offsetof(a0FirmWorkforceDeltaMessage, firmGlobalId);
    offsets[1] = (MPI_Aint)offsetof(a0FirmWorkforceDeltaMessage, delta);

    MPI_Datatype types[2] = { MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG };

    MPI_Type_create_struct(2, blockLengths, offsets, types, &t);
    MPI_Type_commit(&t);
    return t;
}

void BailInBailOut::a0ComputeFirmWorkforceCost(
    unsigned int firmCountTotal,
    unsigned int mpiRank,
    unsigned int mpiSize,
    unsigned int firmGlobalStartIndex,
    unsigned int firmCountForRank,
    const vector<unsigned short>& localWorkerFirmID,
    unsigned int wage,
    vector<uint64_t>& localFirmWorkforceCost
) {
    (void)mpiRank;

    // Local output
    localFirmWorkforceCost.assign(firmCountForRank, 0ULL);

    // Aggregate locally first so we don't have to send a message for
    // every single worker individually
    unordered_map<unsigned int, uint64_t> localAggregator;
    localAggregator.reserve(localWorkerFirmID.size());

    for (unsigned int i = 0; i < localWorkerFirmID.size(); i++) {
        unsigned int firmGlobalId = (unsigned int)localWorkerFirmID[i];
        localAggregator[firmGlobalId] += (uint64_t)wage;
    }

    vector<vector<a0FirmWorkforceDeltaMessage>> toRank(mpiSize);

    for (auto& keyValue : localAggregator) {
        unsigned int firmGlobalId = keyValue.first;
        uint64_t delta = keyValue.second;
        if (delta == 0) {
            continue;
        }

        unsigned int ownerRank = ownerRankFromGlobalId(
            firmCountTotal,
            mpiSize,
            firmGlobalId
        );

        toRank[ownerRank].push_back(a0FirmWorkforceDeltaMessage{
                firmGlobalId,
                delta
            }
        );
    }

    static MPI_Datatype msgType = a0MakeFirmWorkforceDeltaType();

    vector<a0FirmWorkforceDeltaMessage> received;
    allToAllvExchange<a0FirmWorkforceDeltaMessage>(
        mpiSize,
        toRank,
        received,
        msgType
    );

    // Apply received deltas to firms owned by this rank
    for (unsigned int i = 0; i < received.size(); i++) {
        unsigned int firmGlobalId = received[i].firmGlobalId;
        uint64_t delta = received[i].delta;

        unsigned int localFirm = firmGlobalId - firmGlobalStartIndex;
        if (localFirm < firmCountForRank) {
            localFirmWorkforceCost[localFirm] += delta;
        }
        else {
            // If entered,there's an error in mapping 
            printf("A0 Error: received firmGlobalId=%u not owned by this rank\n", firmGlobalId);
        }
    }
}

void BailInBailOut::a6RequestFirmLoans(
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
) {
    for (unsigned int f = 0; f < localFirmLiquidity.size(); f++) {

        // Only request if liquidity is negative
        if (localFirmLiquidity[f] >= 0) {
            continue;
        }

        unsigned int firmGlobalId = firmGlobalStartIndex + f;

        // Request exactly deficit amount
        uint64_t amountRequested = (uint64_t)(-localFirmLiquidity[f]);

        vector<unsigned int>& candidates = localFirmNeighbors[f];

        // If there are no neighbors, it's impossible to request loans
        if (candidates.empty()) {
            continue;
        }

        // Start with first bank to guarantee some value
        unsigned int bestBankGlobalId = candidates[0];

        // Rebuild interest rates from seed to avoid MPI calls
        unsigned short bestRate = getInterestRate(
            baseSeed,
            run,
            bestBankGlobalId,
            minInterestRate,
            maxInterestRate
        );

        for (unsigned int i = 1; i < candidates.size(); i++) {
            unsigned int bankGlobalId = candidates[i];
            // Get the rate for that bank
            // NOTE: We can use seeding to do that since it is deterministic
            // baseSeed and time step are constants for this, and so is
            // the run in relation to the current action
            // and global bank IDs are the ones saved as neighbors, so it's
            // easier to regenerate the rate rather than messaging for it
            unsigned short rate = getInterestRate(
                baseSeed,
                run,
                bankGlobalId,
                minInterestRate,
                maxInterestRate
            );

            if (rate < bestRate ||
                (rate == bestRate && bankGlobalId < bestBankGlobalId)) {
                bestRate = rate;
                bestBankGlobalId = bankGlobalId;
            }
        }

        unsigned int ownerRank = ownerRankFromGlobalId(
            bankCountTotal,
            mpiSize,
            bestBankGlobalId
        );

        firmLoanRequestsToRank[ownerRank].push_back(
            FirmLoanRequest{
                firmGlobalId,
                bestBankGlobalId,
                amountRequested,
                bestRate
            }
        );
    }
}

