#include "pm2_bailinbailout.h"

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

inline void BailInBailOut::buildFirmNeighbors(
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
    unsigned int degree = std::min(firmLenderDegree, bankCountTotal);

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
inline void BailInBailOut::repayFirmLoansProRata(
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

        uint64_t seed = makeSeed(
            baseSeed,
            run,
            0,
            bestBankGlobalId,
            STREAM_INTEREST_RATE_SELECTION
        );

        unsigned short bestRate = (unsigned short)randomInRangeFromSeed
        (seed, minInterestRate, maxInterestRate);

        for (unsigned int i = 1; i < candidates.size(); i++) {
            unsigned int bankGlobalId = candidates[i];

            uint64_t candidateSeed = makeSeed(
                baseSeed,
                run,
                0,
                bankGlobalId,
                STREAM_INTEREST_RATE_SELECTION
            );

            // Get the rate for that bank
            // NOTE: We can use seeding to do that since it is deterministic
            // baseSeed and time step are constants for this, and so is
            // the run in relation to the current action
            // and global bank IDs are the ones saved as neighbors, so it's
            // easier to regenerate the rate rather than messaging for it
            unsigned short rate = (unsigned short)
                randomInRangeFromSeed(candidateSeed, minInterestRate, maxInterestRate);

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

void BailInBailOut::d2ProcessFirmLoanRequests(
    unsigned int firmCountTotal,
    unsigned int bankCountForRank,
    unsigned int bankGlobalStartIndex,
    unsigned int mpiSize,
    unsigned short maxFirmLoanPercent,

    vector<int64_t>& localBankLiquidity,

    vector<FirmLoanRequest>& receivedLoanRequests,

    // Output buffers
    vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank, 
    vector<uint64_t>& bankFirmLoanOutflow 
) {
    // Reset outputs for this step
    for (unsigned int r = 0; r < mpiSize; r++) {
        firmLoanAcceptancesToRank[r].clear();
    }
    fill(bankFirmLoanOutflow.begin(), bankFirmLoanOutflow.end(), 0);

    // Bucket requests by local bank index
    vector<vector<FirmLoanRequest>> requestsByLocalBank(bankCountForRank);

    for (unsigned int i = 0; i < receivedLoanRequests.size(); i++) {
        FirmLoanRequest& req = receivedLoanRequests[i];

        // This rank should own the lender bank
        unsigned int localBankIndex = req.lenderBankGlobalId - bankGlobalStartIndex;

        if (localBankIndex < bankCountForRank) {
            requestsByLocalBank[localBankIndex].push_back(req);
        }

        // Otherwise, ignore. Shouldn't happen assuming D1 is correct... :)
    }

    // Process each local bank's requests 
    for (unsigned int localBank = 0; localBank < bankCountForRank; localBank++) {

        vector<FirmLoanRequest>& reqs = requestsByLocalBank[localBank];
        // Skip if no requests
        if (reqs.empty()) {
            continue;
        }

        // Sort by firms' global IDs
        sort(reqs.begin(), reqs.end(),
            [](const FirmLoanRequest& a, const FirmLoanRequest& b) {
                if (a.firmGlobalId != b.firmGlobalId) {
                    return a.firmGlobalId < b.firmGlobalId;
                }
                if (a.amountRequested != b.amountRequested) {
                    return a.amountRequested < b.amountRequested;
                }
                return a.lenderBankGlobalId < b.lenderBankGlobalId;
            }
        );

        int64_t liqSigned = localBankLiquidity[localBank];
        if (liqSigned <= 0) continue;

        uint64_t liquidity = (uint64_t)liqSigned;

        // Lending cap as a fraction of current liquidity
        uint64_t stepCap = (liquidity * (uint64_t)maxFirmLoanPercent) / 100ULL;

        // Also cannot lend more than available liquidity
        uint64_t remainingToLend = (stepCap < liquidity) ? stepCap : liquidity;

        if (remainingToLend == 0) {
            continue;
        }

        for (unsigned int k = 0; k < reqs.size(); k++) {
            FirmLoanRequest& req = reqs[k];

            if (remainingToLend == 0) {
                break;
            }

            uint64_t grant = req.amountRequested;
            if (grant > remainingToLend) grant = remainingToLend;

            if (grant == 0) continue;

            // Record acceptance to firm owner
            unsigned int firmOwner = ownerRankFromGlobalId(
                firmCountTotal,
                mpiSize,
                req.firmGlobalId
            );

            firmLoanAcceptancesToRank[firmOwner].push_back(
                d2FirmLoanAcceptance{
                    req.firmGlobalId,
                    req.lenderBankGlobalId,
                    grant
                }
            );

            // Track how much this bank is lending this step (apply in D.3)
            bankFirmLoanOutflow[localBank] += grant;

            remainingToLend -= grant;
        }
    }

    // This inbox is one step use
    receivedLoanRequests.clear();
}