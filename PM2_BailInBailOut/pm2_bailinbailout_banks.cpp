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
    MPI_Datatype t;
    int blockLengths[2] = { 1, 1 };
    MPI_Aint offsets[2];
    offsets[0] = (MPI_Aint)offsetof(c1BankDeltaMessage, bankGlobalId);
    offsets[1] = (MPI_Aint)offsetof(c1BankDeltaMessage, delta);
    MPI_Datatype types[2] = { MPI_UNSIGNED, MPI_LONG_LONG };

    MPI_Type_create_struct(2, blockLengths, offsets, types, &t);
    MPI_Type_commit(&t);
    return t;
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
    vector<uint64_t>& bankFirmLoanOutflow,

    // Debug output 
    vector<DDebugGrantEdge>* grantsThisStep,
    // Debug flag
    bool debug
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

            if (debug) {
                grantsThisStep->push_back(
                    DDebugGrantEdge{
                        req.firmGlobalId,
                        req.lenderBankGlobalId,
                        grant
                    }
                );
            }

            // Track how much this bank is lending this step (apply in D.3)
            bankFirmLoanOutflow[localBank] += grant;

            remainingToLend -= grant;
        }
    }

    // This inbox is one step use
    receivedLoanRequests.clear();
}
    
MPI_Datatype BailInBailOut::d4MakeFirmLoanAcceptanceType() {
    MPI_Datatype t;

    int blockLengths[3] = { 1, 1, 1 };
    MPI_Aint offsets[3];
    offsets[0] = (MPI_Aint)offsetof(d2FirmLoanAcceptance, firmGlobalId);
    offsets[1] = (MPI_Aint)offsetof(d2FirmLoanAcceptance, lenderBankGlobalId);
    offsets[2] = (MPI_Aint)offsetof(d2FirmLoanAcceptance, amountGranted);

    MPI_Datatype types[3] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG };

    MPI_Type_create_struct(3, blockLengths, offsets, types, &t);
    MPI_Type_commit(&t);
    return t;
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
    MPI_Datatype t;

    int blockLengths[2] = { 1, 1 };
    MPI_Aint offsets[2];
    offsets[0] = (MPI_Aint)offsetof(e1InterbankRepaymentMessage, lenderBankGlobalId);
    offsets[1] = (MPI_Aint)offsetof(e1InterbankRepaymentMessage, amount);

    MPI_Datatype types[2] = { MPI_UNSIGNED, MPI_UNSIGNED_LONG_LONG };

    MPI_Type_create_struct(2, blockLengths, offsets, types, &t);
    MPI_Type_commit(&t);
    return t;
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

            // Don’t repay more cash than borrower has
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
            cerr << "E1 ERROR: Received lenderGlobalId=%u not owned by this rank\n", lenderGlobalId;
        }
    }
}

void BailInBailOut::e2ApplyInterestOnAllLoans(
    vector<vector<DebtEntry>>& eLocalBankDebts,

) {

   /*
   important to have - 
   debts
   information to call interest rate function
   */
}