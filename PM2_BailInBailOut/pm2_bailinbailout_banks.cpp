#include "pm2_bailinbailout.h"

/*<===========================================================================>

    Winter 2026
    Author: Kyryll Kotyk
    Advisor: Munehiro Fukuda

    This file contains all Bank related logic, including custom MPI datatypes
    and messaging.
    Implementation ownership considers Banks to be dominant.
    When Firms interact with Banks, the implementation is in the Banks file.

    Style Note:
        One parameter per line for function signatures
        Parantheses are closed on a new line, followed by body brace opening
        This rule does not apply if there are no parameters

<===========================================================================>*/

MPI_Datatype BailInBailOut::c1MakeBankDeltaMessageType() {
    MPI_Datatype raw;
    int blockLengths[2] = { 1, 1 };
    MPI_Aint offsets[2];
    offsets[0] = (MPI_Aint)offsetof(c1BankDeltaMessage, bankGlobalId);
    offsets[1] = (MPI_Aint)offsetof(c1BankDeltaMessage, delta);
    MPI_Datatype types[2] = { MPI_UNSIGNED, MPI_LONG_LONG };

    MPI_Type_create_struct(2, blockLengths, offsets, types, &raw);

    MPI_Datatype resized;
    MPI_Type_create_resized(raw, 0, (MPI_Aint)sizeof(c1BankDeltaMessage), &resized);

    MPI_Type_commit(&resized);
    MPI_Type_free(&raw);
    return resized;
}

void BailInBailOut::c1SendMoneyToBanks(
    unsigned int bankCountTotal,
    unsigned int mpiRank,
    unsigned int mpiSize,
    unsigned int bankGlobalStartIndex,
    unsigned int bankCountForRank,
    vector<int64_t>& localBankLiquidity,
    vector<int64_t>& wageDepositBuffer,
    vector<int64_t>& bankIncomingFromFirmRepay
) {
    vector<vector<c1BankDeltaMessage>> toRank(mpiSize);

    for (unsigned int bankGID = 0; bankGID < bankCountTotal; bankGID++) {
        int64_t delta = wageDepositBuffer[bankGID] +
            bankIncomingFromFirmRepay[bankGID];
        if (delta == 0) continue;

        unsigned int ownerRank = ownerRankFromGlobalId(
            bankCountTotal,
            mpiSize,
            bankGID
        );

        toRank[ownerRank].push_back(c1BankDeltaMessage{ bankGID, delta });
    }

    static MPI_Datatype bankDeltaType = c1MakeBankDeltaMessageType();
    vector<c1BankDeltaMessage> recvBuffer;

    allToAllvExchange(mpiSize, toRank, recvBuffer, bankDeltaType);

    for (unsigned int i = 0; i < recvBuffer.size(); i++) {
        unsigned int bankGlobalId = recvBuffer[i].bankGlobalId;
        int64_t delta = recvBuffer[i].delta;

        unsigned int localIndex = bankGlobalId - bankGlobalStartIndex;
        localBankLiquidity[localIndex] += delta;
    }

    fill(
        wageDepositBuffer.begin(),
        wageDepositBuffer.end(),
        0
    );

    fill(
        bankIncomingFromFirmRepay.begin(),
        bankIncomingFromFirmRepay.end(),
        0
    );

}

MPI_Datatype BailInBailOut::d1MakeFirmLoanRequestType() {
    MPI_Datatype raw;

    int blockLengths[4] = { 1,1,1,1 };
    MPI_Aint offsets[4] = {
        (MPI_Aint)offsetof(FirmLoanRequest, firmGlobalId),
        (MPI_Aint)offsetof(FirmLoanRequest, lenderBankGlobalId),
        (MPI_Aint)offsetof(FirmLoanRequest, amountRequested),
        (MPI_Aint)offsetof(FirmLoanRequest, lenderInterestRate),
    };
    MPI_Datatype types[4] = {
        MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_SHORT
    };

    MPI_Type_create_struct(4, blockLengths, offsets, types, &raw);

    MPI_Datatype resized;
    MPI_Type_create_resized(raw, 0, (MPI_Aint)sizeof(FirmLoanRequest), &resized);

    MPI_Type_commit(&resized);
    MPI_Type_free(&raw);

    return resized;
}
/*
MPI_Datatype BailInBailOut::d1MakeFirmLoanRequestType() {
    MPI_Datatype t;

    int blockLengths[4] = { 1, 1, 1, 1 };
    MPI_Aint offsets[4];
    offsets[0] = (MPI_Aint)offsetof(FirmLoanRequest, firmGlobalId);
    offsets[1] = (MPI_Aint)offsetof(FirmLoanRequest, lenderBankGlobalId);
    offsets[2] = (MPI_Aint)offsetof(FirmLoanRequest, amountRequested);
    offsets[3] = (MPI_Aint)offsetof(FirmLoanRequest, lenderInterestRate);

    MPI_Datatype types[4] = {
        MPI_UNSIGNED,
        MPI_UNSIGNED,
        MPI_UNSIGNED_LONG_LONG,
        MPI_UNSIGNED_SHORT
    };

    MPI_Type_create_struct(4, blockLengths, offsets, types, &t);
    MPI_Type_commit(&t);
    return t;
}
*/


void BailInBailOut::d1SendFirmLoanRequests(
    unsigned int mpiSize,
    vector<vector<FirmLoanRequest>>& firmLoanRequestsToRank,
    vector<FirmLoanRequest>& receivedLoanRequests
) {
    // Use the datatype made for this
    static MPI_Datatype firmLoanReqType = d1MakeFirmLoanRequestType();

    // All to all variable sized exchange
    allToAllvExchange(
        mpiSize,
        firmLoanRequestsToRank,
        receivedLoanRequests,
        firmLoanReqType
    );
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

        // Pro-rata allocation: all requesting firms share capacity proportionally
        // instead of first-come-first-served by firm ID
        {
            const size_t n = reqs.size();

            uint64_t totalRequested = 0;
            for (size_t k = 0; k < n; k++) {
                totalRequested += reqs[k].amountRequested;
            }

            vector<uint64_t> grants(n, 0);
            vector<uint64_t> remainders(n, 0);
            uint64_t sumGrants = 0;

            for (size_t k = 0; k < n; k++) {
                unsigned __int128 scaled =
                    (unsigned __int128)remainingToLend * reqs[k].amountRequested;
                grants[k] = (uint64_t)(scaled / totalRequested);
                remainders[k] = (uint64_t)(scaled % totalRequested);
                // Never grant more than requested
                if (grants[k] > reqs[k].amountRequested) {
                    grants[k] = reqs[k].amountRequested;
                }
                sumGrants += grants[k];
            }

            // Distribute leftover 1-by-1 to largest-remainder firms first,
            // breaking ties by lowest firm global ID for determinism
            uint64_t leftover = remainingToLend - sumGrants;
            if (leftover > 0) {
                vector<size_t> order;
                order.reserve(n);
                for (size_t k = 0; k < n; k++) {
                    if (reqs[k].amountRequested > grants[k]) {
                        order.push_back(k);
                    }
                }
                sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                    if (remainders[a] != remainders[b]) {
                        return remainders[a] > remainders[b];
                    }
                    return reqs[a].firmGlobalId < reqs[b].firmGlobalId;
                    });
                for (size_t i = 0; i < order.size() && leftover > 0; i++) {
                    grants[order[i]]++;
                    leftover--;
                }
            }

            // Send acceptances
            for (size_t k = 0; k < n; k++) {
                uint64_t grant = grants[k];
                if (grant == 0) {
                    continue;
                }
                unsigned int firmOwner = ownerRankFromGlobalId(
                    firmCountTotal,
                    mpiSize,
                    reqs[k].firmGlobalId
                );
                firmLoanAcceptancesToRank[firmOwner].push_back(
                    d2FirmLoanAcceptance{
                        reqs[k].firmGlobalId,
                        reqs[k].lenderBankGlobalId,
                        grant
                    }
                );
                bankFirmLoanOutflow[localBank] += grant;
            }
        }
    }

    // This inbox is one step use
    receivedLoanRequests.clear();
}

MPI_Datatype BailInBailOut::d4MakeFirmLoanAcceptanceType() {
    MPI_Datatype raw;

    int blockLengths[3] = { 1, 1, 1 };
    MPI_Aint offsets[3];
    offsets[0] = (MPI_Aint)offsetof(d2FirmLoanAcceptance, firmGlobalId);
    offsets[1] = (MPI_Aint)offsetof(d2FirmLoanAcceptance, lenderBankGlobalId);
    offsets[2] = (MPI_Aint)offsetof(d2FirmLoanAcceptance, amountGranted);

    MPI_Datatype types[3] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG };

    MPI_Type_create_struct(3, blockLengths, offsets, types, &raw);

    MPI_Datatype resized;
    MPI_Type_create_resized(raw, 0, (MPI_Aint)sizeof(d2FirmLoanAcceptance), &resized);

    MPI_Type_commit(&resized);
    MPI_Type_free(&raw);
    return resized;
}

void BailInBailOut::d4SendFirmLoanAcceptances(
    unsigned int mpiSize,
    vector<vector<d2FirmLoanAcceptance>>& firmLoanAcceptancesToRank,
    vector<d2FirmLoanAcceptance>& receivedAcceptances
) {
    // Use the datatype made for this
    static MPI_Datatype firmLoanAccType = d4MakeFirmLoanAcceptanceType();

    // All to all variable sized exchange
    allToAllvExchange(
        mpiSize,
        firmLoanAcceptancesToRank,
        receivedAcceptances,
        firmLoanAccType
    );
}

void BailInBailOut::d5ApplyFirmLoanAcceptances(
    unsigned int firmGlobalStartIndex,
    vector<int64_t>& localFirmLiquidity,
    vector<vector<DebtEntry>>& localFirmDebts,
    vector<d2FirmLoanAcceptance>& receivedAcceptances
) {
    for (unsigned int i = 0; i < receivedAcceptances.size(); i++) {

        d2FirmLoanAcceptance& acc = receivedAcceptances[i];

        // Convert global firm id to local firm index
        unsigned int localFirmIndex = acc.firmGlobalId - firmGlobalStartIndex;

        // If routing is correct, this should always be in range
        if (localFirmIndex >= localFirmLiquidity.size()) {
            continue;
        }

        // Apply cash inflow to firm
        localFirmLiquidity[localFirmIndex] += (int64_t)acc.amountGranted;

        // Record debt on the borrower side
        localFirmDebts[localFirmIndex].push_back(
            DebtEntry{
                acc.lenderBankGlobalId,
                acc.amountGranted
            }
        );
    }

    // Inbox is one step use
    receivedAcceptances.clear();
}

MPI_Datatype BailInBailOut::e1MakeInterbankRepaymentType() {
    MPI_Datatype raw;

    int blockLengths[2] = { 1, 1 };
    MPI_Aint offsets[2];
    offsets[0] = (MPI_Aint)offsetof(e1InterbankRepaymentMessage, lenderBankGlobalId);
    offsets[1] = (MPI_Aint)offsetof(e1InterbankRepaymentMessage, amount);

    MPI_Datatype types[2] = { MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG };

    MPI_Type_create_struct(2, blockLengths, offsets, types, &raw);

    MPI_Datatype resized;
    MPI_Type_create_resized(raw, 0, (MPI_Aint)sizeof(e1InterbankRepaymentMessage), &resized);

    MPI_Type_commit(&resized);
    MPI_Type_free(&raw);

    return resized;
}

void BailInBailOut::e1RepayInterbankLoans(
    unsigned int bankCountTotal,
    unsigned int mpiRank,
    unsigned int mpiSize,
    unsigned int bankGlobalStartIndex,
    unsigned int bankCountForRank,
    unsigned short bankRepayPercent,
    vector<int64_t>& localBankLiquidity,
    vector<vector<DebtEntry>>& localBankDebts
) {
    (void)mpiRank;

    // Build outgoing credit messages to lenders per destination rank
    vector<vector<e1InterbankRepaymentMessage>> toRank(mpiSize);

    for (unsigned int localBank = 0; localBank < bankCountForRank; localBank++) {

        // Don't repa if negative liquidity
        if (localBankLiquidity[localBank] <= 0) {
            continue;
        }

        vector<DebtEntry>& debts = localBankDebts[localBank];

        // Skip if nothing to repay
        if (debts.empty()) {
            continue;
        }

        for (unsigned int i = 0; i < debts.size(); i++) {
            DebtEntry& entry = debts[i];

            // Skip this debt entry if it's 0
            if (entry.amount == 0) {
                continue;
            }

            // bankRepayPercent% of remaining debt
            uint64_t repay =
                (uint64_t)((entry.amount * (uint64_t)bankRepayPercent) / 100ULL);

            if (repay == 0) {
                continue;
            }

            // Don't repay more cash than borrower has
            uint64_t cash =
                (localBankLiquidity[localBank] > 0) ?
                (uint64_t)localBankLiquidity[localBank] : 0ULL;

            // If there's not enough to pay off the desired amount,
            // pay all of the available cash
            if (repay > cash) {
                repay = cash;
            }

            // If there's nothing left to repay, stop
            if (repay == 0) {
                break;
            }

            // Borrower subtracts instantly 
            localBankLiquidity[localBank] -= (int64_t)repay;

            // Reduce the borrower's stored liability
            entry.amount -= (uint64_t)repay;

            // Send credit to lender owner rank similarly to C.1
            unsigned int lenderGlobalId = entry.lenderBankGlobalId;

            unsigned int ownerRank = ownerRankFromGlobalId(
                bankCountTotal,
                mpiSize,
                lenderGlobalId
            );

            toRank[ownerRank].push_back(e1InterbankRepaymentMessage{
                lenderGlobalId,
                repay
                });

            // If borrower cash is now below 0, stop repaying further debts
            if (localBankLiquidity[localBank] <= 0) {
                break;
            }
        }

        // Remove fully repaid entries to save space
        unsigned int write = 0;
        for (unsigned int i = 0; i < debts.size(); i++) {
            if (debts[i].amount != 0) {
                debts[write++] = debts[i];
            }
        }
        debts.resize(write);
    }

    // Exchange and apply credits on lender owners 
    static MPI_Datatype msgType = e1MakeInterbankRepaymentType();

    vector<e1InterbankRepaymentMessage> received;
    allToAllvExchange<e1InterbankRepaymentMessage>(
        mpiSize,
        toRank,
        received,
        msgType
    );

    // Apply lender credits to banks owned by this rank
    for (unsigned int i = 0; i < received.size(); i++) {
        unsigned int lenderGlobalId = received[i].lenderBankGlobalId;
        uint64_t amount = received[i].amount;

        unsigned int lenderLocal = lenderGlobalId - bankGlobalStartIndex;
        if (lenderLocal < bankCountForRank) {
            localBankLiquidity[lenderLocal] += (int64_t)amount;
        }
        else {
            // Should NEVER be entered. Indicates mismatch in range mapping.
            printf("E1 ERROR: Received lenderGlobalId=%u not owned by this rank\n", lenderGlobalId);
        }
    }
}

void BailInBailOut::e2ApplyInterestOnAllLoans(
    vector<vector<DebtEntry>>& eLocalBankDebts,
    vector<vector<DebtEntry>>& localFirmDebts,
    uint64_t baseSeed,
    unsigned int run,
    uint64_t minVal,
    uint64_t maxVal
) {

    // Go through every BANK debt entry
    for (int i = 0; i < eLocalBankDebts.size(); i++) {
        for (DebtEntry& entry : eLocalBankDebts[i]) {
            entry.amount = (entry.amount * (100 +
                getInterestRate(
                    baseSeed,
                    run,
                    entry.lenderBankGlobalId,
                    minVal,
                    maxVal
                ))) / 100ULL;
        }
    }

    // Go through every FIRM debt entry
    for (int i = 0; i < localFirmDebts.size(); i++) {
        for (DebtEntry& entry : localFirmDebts[i]) {
            entry.amount = (entry.amount * (100 +
                getInterestRate(
                    baseSeed,
                    run,
                    entry.lenderBankGlobalId,
                    minVal,
                    maxVal
                ))) / 100ULL;
        }
    }

}

MPI_Datatype BailInBailOut::e3MakeInterbankLoanRequestType() {
    MPI_Datatype raw;

    int blockLengths[4] = { 1, 1, 1, 1 };
    MPI_Aint offsets[4] = {
        (MPI_Aint)offsetof(e3InterbankLoanRequest, borrowerBankGlobalId),
        (MPI_Aint)offsetof(e3InterbankLoanRequest, lenderBankGlobalId),
        (MPI_Aint)offsetof(e3InterbankLoanRequest, amountRequested),
        (MPI_Aint)offsetof(e3InterbankLoanRequest, lenderInterestRate),
    };
    MPI_Datatype types[4] = {
        MPI_UNSIGNED,
        MPI_UNSIGNED,
        MPI_UNSIGNED_LONG_LONG,
        MPI_UNSIGNED_SHORT
    };

    MPI_Type_create_struct(4, blockLengths, offsets, types, &raw);

    MPI_Datatype resized;
    MPI_Type_create_resized(raw, 0, (MPI_Aint)sizeof(e3InterbankLoanRequest), &resized);

    MPI_Type_commit(&resized);
    MPI_Type_free(&raw);

    return resized;
}

MPI_Datatype BailInBailOut::e3MakeInterbankLoanAcceptanceType() {
    MPI_Datatype raw;

    int blockLengths[3] = { 1, 1, 1 };
    MPI_Aint offsets[3] = {
        (MPI_Aint)offsetof(e3InterbankLoanAcceptance, borrowerBankGlobalId),
        (MPI_Aint)offsetof(e3InterbankLoanAcceptance, lenderBankGlobalId),
        (MPI_Aint)offsetof(e3InterbankLoanAcceptance, amountGranted)
    };
    MPI_Datatype types[3] = {
        MPI_UNSIGNED,
        MPI_UNSIGNED,
        MPI_UNSIGNED_LONG_LONG
    };

    MPI_Type_create_struct(3, blockLengths, offsets, types, &raw);

    MPI_Datatype resized;
    MPI_Type_create_resized(raw, 0, (MPI_Aint)sizeof(e3InterbankLoanAcceptance), &resized);

    MPI_Type_commit(&resized);
    MPI_Type_free(&raw);

    return resized;
}

MPI_Datatype BailInBailOut::e3MakeNeighborLiquidityMessageType() {
    MPI_Datatype raw;

    int blockLengths[2] = { 1, 1 };
    MPI_Aint offsets[2] = {
        (MPI_Aint)offsetof(e3NeighborLiquidityMessage, bankGlobalId),
        (MPI_Aint)offsetof(e3NeighborLiquidityMessage, liquidity)
    };
    MPI_Datatype types[2] = {
        MPI_UNSIGNED,
        MPI_LONG_LONG
    };

    MPI_Type_create_struct(2, blockLengths, offsets, types, &raw);

    MPI_Datatype resized;
    MPI_Type_create_resized(raw, 0, (MPI_Aint)sizeof(e3NeighborLiquidityMessage), &resized);

    MPI_Type_commit(&resized);
    MPI_Type_free(&raw);

    return resized;
}

void BailInBailOut::e3BorrowInterbankIfNegativeLiquidity(
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
) {
    // Due to the length of this function, it follows a subphase format

    // <======================================================================>
    // <==================E.3.1: Neighbor liquidity exchange==================>
    // <======================================================================>

    // Build messages to send
    // For each local bank, send its liquidity to all neighbors
    // Route messages by owner rank of the neighbor, and send only to ranks that
    // need this information

    vector<vector<e3NeighborLiquidityMessage>> liquidityMessagesToRank(mpiSize);

    // Helper will compute whether each rank will send anything based on
    // whether liquidityMessagesToRank[rank] is empty or not

    // Iterate over each bank on this rank
    for (unsigned int localBank = 0; localBank < bankCountForRank; localBank++) {

        // Find the local bank's global ID
        unsigned int bankGlobalId = bankGlobalStartIndex + localBank;
        // Get its liquidity
        int64_t bankLiquidity = localBankLiquidity[localBank];

        // Liquidity to each neighbor
        vector<unsigned int>& neighbors = localBankNeighbors[localBank];

        // Go through all the neighbors
        for (unsigned int i = 0; i < neighbors.size(); i++) {
            // Get the neighbor's global ID
            unsigned int neighborGlobalId = neighbors[i];

            // Get the specific rank it's on
            unsigned int neighborOwnerRank = ownerRankFromGlobalId(
                bankCountTotal,
                mpiSize,
                neighborGlobalId
            );

            // If neighbor is local, no message needed
            if (neighborOwnerRank == mpiRank) {
                continue;
            }

            // Add the message to the accumulator
            liquidityMessagesToRank[neighborOwnerRank].push_back(
                e3NeighborLiquidityMessage{
                    bankGlobalId,
                    bankLiquidity
                }
            );
        }
    }

    static MPI_Datatype liquidityType = e3MakeNeighborLiquidityMessageType();

    vector<e3NeighborLiquidityMessage> receivedNeighborLiquidity;

    // Direct sparse exchange for neighbor liquidity
    sparseDirectExchange<e3NeighborLiquidityMessage>(
        mpiRank,
        mpiSize,
        liquidityMessagesToRank,
        receivedNeighborLiquidity,
        liquidityType,
        /*countTag*/ 3201,
        /*dataTag*/  3202
    );

    // Build a liquidity view for banks that might be considered as lenders
    // Local banks from localBankLiquidity
    // Remote banks from neighbor liquidity messages
    vector<int64_t> lenderLiquidityView(bankCountTotal, INT64_MIN);

    auto pushLocalLiquidityIntoView = [&]() {
        for (unsigned int localBank = 0; localBank < bankCountForRank; localBank++) {
            unsigned int globalId = bankGlobalStartIndex + localBank;
            lenderLiquidityView[globalId] = localBankLiquidity[localBank];
        }
        };

    pushLocalLiquidityIntoView();

    // Fill remote liquidity from incoming messages
    for (unsigned int i = 0; i < receivedNeighborLiquidity.size(); i++) {
        lenderLiquidityView[receivedNeighborLiquidity[i].bankGlobalId] =
            receivedNeighborLiquidity[i].liquidity;
    }

    // <======================================================================>
    // <=================E.3.2: Build interbank loan requests=================>
    // <======================================================================>

    vector<vector<e3InterbankLoanRequest>> requestsToRank(mpiSize);

    auto getLenderCapFromLiquidity = [&](int64_t lenderLiquidity) -> uint64_t {
        if (lenderLiquidity <= 0) {
            return 0ULL;
        }
        uint64_t liq = (uint64_t)lenderLiquidity;
        return (liq * (uint64_t)maxInterbankLoanPercent) / 100ULL;
        };

    for (unsigned int localBorrower = 0; localBorrower < bankCountForRank; localBorrower++) {
        // If liquidity is not negative, don't need loans
        if (localBankLiquidity[localBorrower] >= 0) {
            continue;
        }

        // Get the global ID of the borrower bank
        unsigned int borrowerGlobalId = bankGlobalStartIndex + localBorrower;

        // Figure out the deficit
        // Deficit is how much is needed to be loaned
        uint64_t deficitRemaining = (uint64_t)(-localBankLiquidity[localBorrower]);

        // Extra check just to make sure there's some deficit
        if (deficitRemaining == 0) {
            continue;
        }

        // All the candidates for loaner
        vector<unsigned int>& candidates = localBankNeighbors[localBorrower];

        // If this bank how no loaner candidates, skip
        if (candidates.empty()) {
            continue;
        }

        unsigned int sampleCount = (unsigned int)maxInterbankLenderSamplingK;

        // If too little nodes to reach the max samplers, then
        // set the sampling size to total candidate count
        if (sampleCount > (unsigned int)candidates.size()) {
            sampleCount = (unsigned int)candidates.size();
        }

        vector<e3LenderOption> lenderOptions;
        lenderOptions.reserve(sampleCount);

        for (unsigned int i = 0; i < sampleCount; i++) {
            unsigned int lenderGlobalId = candidates[i];

            // Interest rate of this specific lender
            unsigned short rate = getInterestRate(
                baseSeed,
                run,
                lenderGlobalId,
                minInterestRate,
                maxInterestRate
            );

            // Liquidity of this specific lender (INT64_MIN means unknown)
            int64_t liquidity = lenderLiquidityView[lenderGlobalId];

            lenderOptions.push_back(e3LenderOption{
                lenderGlobalId,
                rate,
                liquidity
                });
        }

        // Sort lenders by:
        // Lowest interest rate
        // If tied, higher liquidity
        // If STILL tied, lowest global ID
        sort(lenderOptions.begin(), lenderOptions.end(),
            [](const e3LenderOption& a, const e3LenderOption& b) {
                if (a.interestRate != b.interestRate) {
                    return a.interestRate < b.interestRate;
                }
                if (a.liquidity != b.liquidity) {
                    return a.liquidity > b.liquidity;
                }
                return a.lenderGlobalId < b.lenderGlobalId;
            }
        );

        // Request loans from multiple lenders if needed
        // Borrow as much as possible from best lender
        // If not enough, go to next, repeat
        for (unsigned int i = 0; i < lenderOptions.size(); i++) {
            if (deficitRemaining == 0) {
                break;
            }

            unsigned int lenderGlobalId = lenderOptions[i].lenderGlobalId;

            // If lender liquidity is unknown, still possible to request,
            // but it's hard to cap properly
            // In that case, request the remaining deficit, and let lender
            // handle whether it can grant that request or not
            int64_t lenderLiquidity = lenderOptions[i].liquidity;

            uint64_t plannedRequest = deficitRemaining;

            // Make sure the lender has some liquidity
            if (lenderLiquidity > 0) {
                // Planned request is capped by estimated lender capacity
                uint64_t estimatedCap = getLenderCapFromLiquidity(lenderLiquidity);

                if (plannedRequest > estimatedCap) {
                    plannedRequest = estimatedCap;
                }
            }

            // If estimated cap is 0 (or lender had negative liquidity),
            // skip requesting this lender
            if (plannedRequest == 0) {
                continue;
            }

            // Get the owner rank of the lender
            unsigned int lenderOwnerRank = ownerRankFromGlobalId(
                bankCountTotal,
                mpiSize,
                lenderGlobalId
            );

            // Add this to pool of requests for this rank
            requestsToRank[lenderOwnerRank].push_back(e3InterbankLoanRequest{
                borrowerGlobalId,
                lenderGlobalId,
                plannedRequest,
                lenderOptions[i].interestRate
                });

            // Assume attempt to cover this portion of the deficit using this lender
            // If lender cannot fully grant, borrower will remain negative and
            // try again next timestep
            deficitRemaining -= plannedRequest;
        }
    }

    // <======================================================================>
    // <==============E.3.3: Direct sparse exchange for requests==============>
    // <======================================================================>

    static MPI_Datatype reqType = e3MakeInterbankLoanRequestType();

    vector<e3InterbankLoanRequest> receivedRequests;

    // Direct sparse exchange for requests
    sparseDirectExchange<e3InterbankLoanRequest>(
        mpiRank,
        mpiSize,
        requestsToRank,
        receivedRequests,
        reqType,
        /*countTag*/ 3101,
        /*dataTag*/  3102
    );

    // ------------------------------------------------------------------
    // Lenders process requests deterministically and grant up to capacity
    // ------------------------------------------------------------------

    vector<vector<e3InterbankLoanRequest>> reqsByLocalLender(bankCountForRank);

    for (unsigned int i = 0; i < receivedRequests.size(); i++) {
        e3InterbankLoanRequest& req = receivedRequests[i];

        unsigned int localLender = req.lenderBankGlobalId - bankGlobalStartIndex;
        if (localLender < bankCountForRank) {
            reqsByLocalLender[localLender].push_back(req);
        }
    }

    vector<vector<e3InterbankLoanAcceptance>> acceptancesToRank(mpiSize);

    for (unsigned int localLender = 0; localLender < bankCountForRank; localLender++) {

        vector<e3InterbankLoanRequest>& reqs = reqsByLocalLender[localLender];
        if (reqs.empty()) {
            continue;
        }

        // Sort by borrower global ID for determinism
        sort(reqs.begin(), reqs.end(),
            [](const e3InterbankLoanRequest& a, const e3InterbankLoanRequest& b) {
                if (a.borrowerBankGlobalId != b.borrowerBankGlobalId) {
                    return a.borrowerBankGlobalId < b.borrowerBankGlobalId;
                }
                if (a.amountRequested != b.amountRequested) {
                    return a.amountRequested < b.amountRequested;
                }
                return a.lenderBankGlobalId < b.lenderBankGlobalId;
            }
        );

        int64_t lenderLiquiditySigned = localBankLiquidity[localLender];
        if (lenderLiquiditySigned <= 0) {
            continue;
        }

        uint64_t lenderLiquidity = (uint64_t)lenderLiquiditySigned;

        // Max lending capacity this step
        uint64_t capFrac = (lenderLiquidity * (uint64_t)maxInterbankLoanPercent) / 100ULL;
        uint64_t remainingToLend = (capFrac < lenderLiquidity) ? capFrac : lenderLiquidity;

        if (remainingToLend == 0) {
            continue;
        }

        // Pro-rata allocation: all requesting borrowers share capacity
        // proportionally instead of first-come-first-served by borrower ID
        {
            const size_t n = reqs.size();

            uint64_t totalRequested = 0;
            for (size_t k = 0; k < n; k++) {
                totalRequested += reqs[k].amountRequested;
            }

            vector<uint64_t> grants(n, 0);
            vector<uint64_t> remainders(n, 0);
            uint64_t sumGrants = 0;

            for (size_t k = 0; k < n; k++) {
                unsigned __int128 scaled =
                    (unsigned __int128)remainingToLend * reqs[k].amountRequested;
                grants[k] = (uint64_t)(scaled / totalRequested);
                remainders[k] = (uint64_t)(scaled % totalRequested);
                // Never grant more than requested
                if (grants[k] > reqs[k].amountRequested) {
                    grants[k] = reqs[k].amountRequested;
                }
                sumGrants += grants[k];
            }

            // Distribute leftover 1-by-1 to largest-remainder borrowers first,
            // breaking ties by lowest borrower global ID for determinism
            uint64_t leftover = remainingToLend - sumGrants;
            if (leftover > 0) {
                vector<size_t> order;
                order.reserve(n);
                for (size_t k = 0; k < n; k++) {
                    if (reqs[k].amountRequested > grants[k]) {
                        order.push_back(k);
                    }
                }
                sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                    if (remainders[a] != remainders[b]) {
                        return remainders[a] > remainders[b];
                    }
                    return reqs[a].borrowerBankGlobalId < reqs[b].borrowerBankGlobalId;
                    });
                for (size_t i = 0; i < order.size() && leftover > 0; i++) {
                    grants[order[i]]++;
                    leftover--;
                }
            }

            // Apply grants: lender subtracts, send acceptances to borrowers
            for (size_t k = 0; k < n; k++) {
                uint64_t grant = grants[k];
                if (grant == 0) {
                    continue;
                }
                localBankLiquidity[localLender] -= (int64_t)grant;

                unsigned int borrowerOwnerRank = ownerRankFromGlobalId(
                    bankCountTotal,
                    mpiSize,
                    reqs[k].borrowerBankGlobalId
                );
                acceptancesToRank[borrowerOwnerRank].push_back(e3InterbankLoanAcceptance{
                    reqs[k].borrowerBankGlobalId,
                    reqs[k].lenderBankGlobalId,
                    grant
                    });
            }
        }
    }

    // ------------------------------------------------------------------
    // Direct sparse exchange for acceptances (same as requests)
    // ------------------------------------------------------------------

    static MPI_Datatype accType = e3MakeInterbankLoanAcceptanceType();

    vector<e3InterbankLoanAcceptance> receivedAcceptances;

    // Direct sparse exchange for acceptances
    sparseDirectExchange<e3InterbankLoanAcceptance>(
        mpiRank,
        mpiSize,
        acceptancesToRank,
        receivedAcceptances,
        accType,
        /*countTag*/ 3111,
        /*dataTag*/  3112
    );

    // ------------------------------------------------------------------
    // Borrowers apply granted loans (cash inflow + record debt)
    // ------------------------------------------------------------------

    for (unsigned int i = 0; i < receivedAcceptances.size(); i++) {
        e3InterbankLoanAcceptance& acc = receivedAcceptances[i];

        unsigned int localBorrower = acc.borrowerBankGlobalId - bankGlobalStartIndex;
        if (localBorrower >= bankCountForRank) {
            continue;
        }

        localBankLiquidity[localBorrower] += (int64_t)acc.amountGranted;

        localBankDebts[localBorrower].push_back(DebtEntry{
            acc.lenderBankGlobalId,
            acc.amountGranted
            });
    }
}

void BailInBailOut::f1f2UpdateBankDistressAndApplyIntervention(
    unsigned int bankCountForRank,
    unsigned int interventionDelay,
    unsigned int policy,
    unsigned short bailInCoveragePercent,
    unsigned short bailOutCoveragePercent,
    vector<int64_t>& localBankLiquidity,
    vector<vector<DebtEntry>>& localBankDebts,
    vector<unsigned int>& bankDistressCount,
    uint64_t& bailMoneyThisStep,
    uint64_t& graceMoneyThisStep
) {
    bailMoneyThisStep = 0;
    graceMoneyThisStep = 0;


    for (unsigned int localBank = 0; localBank < bankCountForRank; localBank++) {
        vector<DebtEntry>& debts = localBankDebts[localBank];
        // Get the total debts of this bank
        uint64_t totalDebtBefore = sumDebtU64(debts);

        bool isInsolvent = localBankLiquidity[localBank] < totalDebtBefore;

        // F.1: Update distress counter
        if (isInsolvent) {
            bankDistressCount[localBank]++;
        }
        // Reset if it's solvent
        else {
            bankDistressCount[localBank] = 0;
            continue;
        }

        // F.2: Apply intervention every interventionDelay distress timesteps
        unsigned int distressCount = bankDistressCount[localBank];

        bool shouldIntervene =
            (distressCount >= interventionDelay) &&
            (distressCount % interventionDelay == 0);
        // If it's not time to intervene yet, don't proceed
        if (!shouldIntervene) {
            continue;
        }

        uint64_t deficitBefore = totalDebtBefore - localBankLiquidity[localBank];

        // Policy 0 -> Bail-In
        if (policy == 0) {
            // Check if there is some debt
            if (totalDebtBefore > 0) {
                uint64_t targetDebtRemoval =
                    percentFloorU64(deficitBefore, bailInCoveragePercent);

                // If the target removal is greater than what is needed,
                // only cover the total debt 
                if (targetDebtRemoval > totalDebtBefore) {
                    targetDebtRemoval = totalDebtBefore;
                }

                if (targetDebtRemoval > 0) {
                    const size_t debtCount = debts.size();

                    // Amount of reduction applied to each debt entry
                    vector<uint64_t> debtReduction(debtCount, 0);
                    // Remainder used to distribute leftover reductions
                    vector<uint64_t> remainderByDebt(debtCount, 0);
                    // Total reduction assigned during the allocation pass
                    uint64_t assignedReduction = 0;

                    // Compute proportional reduction for each debt based on its share
                    // of total interbank liabilities.
                    for (size_t i = 0; i < debtCount; i++) {

                        // Skip debts that are already zero
                        if (debts[i].amount == 0) {
                            continue;
                        }

                        // Use 128 bit to avoid overflow when multiplying
                        unsigned __int128 scaled =
                            (unsigned __int128)targetDebtRemoval *
                            (unsigned __int128)debts[i].amount;

                        // Determine the proportional reduction
                        uint64_t reduction =
                            (uint64_t)(scaled / totalDebtBefore);
                        // Store remainder for distribution
                        uint64_t remainder =
                            (uint64_t)(scaled % totalDebtBefore);

                        // Never remove more than the actual debt amount
                        if (reduction > debts[i].amount) {
                            reduction = debts[i].amount;
                        }

                        debtReduction[i] = reduction;
                        remainderByDebt[i] = remainder;
                        assignedReduction += reduction;
                    }

                    // Remaining reduction still to be assigned
                    uint64_t leftover = targetDebtRemoval - assignedReduction;

                    if (leftover > 0) {
                        // Order of debts eligible for receiving the leftovers
                        vector<size_t> order;
                        order.reserve(debtCount);
                        // Only debts that still have remaining balance are eligible
                        for (size_t i = 0; i < debtCount; i++) {
                            if (debts[i].amount > debtReduction[i]) {
                                order.push_back(i);
                            }
                        }

                        // Sort by remainder (largest first) 
                        sort(
                            order.begin(),
                            order.end(),
                            [&](size_t a, size_t b) {
                                if (remainderByDebt[a] != remainderByDebt[b]) {
                                    return remainderByDebt[a] > remainderByDebt[b];
                                }
                                if (debts[a].lenderBankGlobalId !=
                                    debts[b].lenderBankGlobalId) {
                                    return debts[a].lenderBankGlobalId <
                                        debts[b].lenderBankGlobalId;
                                }
                                return a < b;
                            }
                        );

                        // Distribute leftover reduction one unit at a time
                        for (size_t k = 0; k < order.size() && leftover > 0; k++) {
                            size_t i = order[k];

                            if (debts[i].amount > debtReduction[i]) {
                                debtReduction[i] += 1;
                                leftover -= 1;
                            }
                        }
                    }

                    // Apply the reductions to the actual debt entries
                    uint64_t actuallyRemoved = 0;

                    for (size_t i = 0; i < debtCount; i++) {
                        uint64_t reduction = debtReduction[i];

                        if (reduction == 0) {
                            continue;
                        }

                        debts[i].amount -= reduction;
                        actuallyRemoved += reduction;
                    }

                    // Remove any debts that were fully eliminated
                    for (size_t i = 0; i < debts.size(); ) {
                        if (debts[i].amount == 0) {
                            debts[i] = debts.back();
                            debts.pop_back();
                        }
                        else {
                            i++;
                        }
                    }

                    // Track total debt removed by bail-in this timestep
                    bailMoneyThisStep += actuallyRemoved;
                }
            }
        }

        // Policy 1 -> Bail-Out
        else {
            if (bailOutCoveragePercent > 0) {
                // Find the injection amount as percentage of deficit
                uint64_t injection =
                    percentFloorU64(deficitBefore, bailOutCoveragePercent);
                // Add that cash
                localBankLiquidity[localBank] += (int64_t)injection;
                // Track the liquidity 
                bailMoneyThisStep += injection;
            }
        }

        // F.5: Grace money if still insolvent after intervention
        uint64_t totalDebtAfter = sumDebtU64(debts);

        bool stillInsolvent =
            ((long double)localBankLiquidity[localBank] < (long double)totalDebtAfter);

        if (stillInsolvent) {
            uint64_t deficitAfter =
                (uint64_t)((long double)totalDebtAfter -
                    (long double)localBankLiquidity[localBank]);

            graceMoneyThisStep += deficitAfter;
        }
    }
}
