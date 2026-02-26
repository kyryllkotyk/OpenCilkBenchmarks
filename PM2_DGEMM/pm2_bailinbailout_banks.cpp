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


struct BailInBailOut::c1BankDeltaMessage {
    unsigned int bankGlobalId;
    int64_t delta;
};

inline MPI_Datatype BailInBailOut::c1MakeBankDeltaMessageType() {
    MPI_Datatype t;
    int blockLengths[2] = { 1, 1 };
    MPI_Aint offsets[2];
    offsets[0] = (MPI_Aint)offsetof(BankDeltaMessage, bankGlobalId);
    offsets[1] = (MPI_Aint)offsetof(BankDeltaMessage, delta);
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
    // Combine the two global ID buffers into sparse messages by destination rank
    vector<vector<c1BankDeltaMessage>> toRank(mpiSize);

    for (unsigned int bankGID = 0; bankGID < bankCountTotal; bankGID++) {
        int64_t delta = wageDepositBuffer[bankGID] +
            bankIncomingFromFirmRepay[bankGID];

        // Skip if there's nothing to send
        if (delta == 0) {
            continue;
        }

        // Get the owner rank from global
        unsigned int ownerRank = ownerRankFromGlobalId(
            bankCountTotal,
            mpiSize, 
            bankGID
        );

        toRank[ownerRank].push_back(
            c1BankDeltaMessage{
                bankGID,
                delta
            }
        );
    }

    // Send counts
    vector<int> sendCounts(mpiSize, 0);

    for (unsigned int r = 0; r < mpiSize; r++) {
        sendCounts[r] = (int)toRank[r].size();
    }

    vector<int> recvCounts(mpiSize, 0);
    // MPI command for every node to send info to every other node
    // Send metadata
    MPI_Alltoall(
        sendCounts.data(),
        1,
        MPI_INT,
        recvCounts.data(),
        1,
        MPI_INT,
        MPI_COMM_WORLD
    );

    // Displacements
    vector<int> sendDisplacements(mpiSize, 0), receiveDisplacements(mpiSize, 0);
    int sendTotal = 0, recvTotal = 0;

    for (unsigned int r = 0; r < mpiSize; r++) {
        sendDisplacements[r] = sendTotal;
        sendTotal += sendCounts[r];
        receiveDisplacements[r] = recvTotal;
        recvTotal += recvCounts[r];
    }

    // Flatten send buffer
    vector<c1BankDeltaMessage> sendBuf;
    sendBuf.reserve(sendTotal);
    for (unsigned int r = 0; r < mpiSize; r++) {
        for (unsigned int i = 0; i < toRank[r].size(); i++) {
            sendBuf.push_back(toRank[r][i]);
        }
    }

    vector<c1BankDeltaMessage> recvBuf(recvTotal);

    static MPI_Datatype bankDeltaType = c1MakeBankDeltaMessageType();

    MPI_Alltoallv(
        sendBuf.data(),
        sendCounts.data(),
        sendDisplacements.data(),
        bankDeltaType,
        recvBuf.data(),
        recvCounts.data(),
        receiveDisplacements.data(),
        bankDeltaType,
        MPI_COMM_WORLD
    );

    // Apply received deltas to local bank liquidity
    for (int i = 0; i < recvTotal; i++) {
        unsigned int bankGlobalId = recvBuf[i].bankGlobalId;
        int64_t delta = recvBuf[i].delta;

        unsigned int localIndex = bankGlobalId - bankGlobalStartIndex;
        localBankLiquidity[localIndex] += delta;
    }

    // Reset buffers
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

inline MPI_Datatype BailInBailOut::d1MakeFirmLoanRequestType() {
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

void BailInBailOut::d1SendFirmLoanRequests(
    unsigned int mpiSize,
    vector<vector<FirmLoanRequest>>& firmLoanRequestsToRank,
    vector<FirmLoanRequest>& receivedLoanRequests
) {

    // Build sendCounts
    vector<int> sendCounts(mpiSize, 0);
    for (unsigned int r = 0; r < mpiSize; r++) {
        sendCounts[r] = (int)firmLoanRequestsToRank[r].size();
    }

    // Exchange counts to size receive buffer
    vector<int> recvCounts(mpiSize, 0);
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Build displacements
    vector<int> sendDisplacements(mpiSize, 0);
    vector<int> receiveDisplacements(mpiSize, 0);

    int sendTotal = 0;
    int recvTotal = 0;

    for (unsigned int r = 0; r < mpiSize; r++) {
        sendDisplacements[r] = sendTotal;
        sendTotal += sendCounts[r];

        receiveDisplacements[r] = recvTotal;
        recvTotal += recvCounts[r];
    }

    // Flatten send buffer in rank order to match sendDisplaycements
    vector<FirmLoanRequest> sendBuf;
    sendBuf.reserve(sendTotal);

    for (unsigned int r = 0; r < mpiSize; r++) {
        for (unsigned int i = 0; i < firmLoanRequestsToRank[r].size(); i++) {
            sendBuf.push_back(firmLoanRequestsToRank[r][i]);
        }
    }

    // Receive buffer
    receivedLoanRequests.clear();
    receivedLoanRequests.resize(recvTotal);

    // Use the datatype made for this
    static MPI_Datatype firmLoanReqType = d1MakeFirmLoanRequestType();

    // All to all variable sized exchange
    MPI_Alltoallv(
        sendBuf.data(),
        sendCounts.data(),
        sendDisplacements.data(),
        firmLoanReqType,
        receivedLoanRequests.data(),
        recvCounts.data(),
        receiveDisplacements.data(),
        firmLoanReqType,
        MPI_COMM_WORLD
    );

    // Clear outgoing buffers for next timestep
    for (unsigned int r = 0; r < mpiSize; r++) {
        firmLoanRequestsToRank[r].clear();
    }
}

