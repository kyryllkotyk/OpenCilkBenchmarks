#include "pm2dgemm.h"

inline pair<int, int> split1D(int totalDimensionSize, int processesInDimension, int index) {
    int base = totalDimensionSize / processesInDimension;
    // First extra blocks get +1 extra  
    int extra = totalDimensionSize % processesInDimension;
    int count = base + (index < extra ? 1 : 0);
    int start = index * base + (index < extra ? index : extra);
    return { start, count };
}

inline double valueAt(int i, int j, int seed) {
    uint64_t x = seed;
    x ^= uint64_t(i) * 0x9e3779b97f4a7c15ULL;
    x ^= uint64_t(j) * 0xbf58476d1ce4e5b9ULL;

    x ^= (x >> 30);
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= (x >> 27);
    x *= 0x94d049bb133111ebULL;
    x ^= (x >> 31);

    return (x >> 11) * (1.0 / (1ULL << 53));
}

vector<vector<double>> PM2_DGEMM::runBenchmark(const short runs,
    int aSize, int bSize,
    uint64_t baseSeed) {

    if (aSize != bSize) {
        cerr << "Cannot multiply, sizes must match";
        return {};
    }

    // Figure out the rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Figure out how many ranks there are
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Turn the process list into a 2D array
    int dims[2] = { 0, 0 };
    MPI_Dims_create(size, 2, dims);

    if (dims[0] != dims[1]) {
        cerr << "MPI created a non-square process grid."
            << " Ensure that the ranks (-np #) is a perfect square.\n";
    }

    int processArrayOrder = dims[0];

    // Figure out grid dimensions
    int myRow = rank / processArrayOrder;
    int myCol = rank % processArrayOrder;

    // Find owned dimensions
    pair<int, int> aLocalRows = split1D(aSize, processArrayOrder, myRow);
    pair<int, int> aLocalCols = split1D(aSize, processArrayOrder, myCol);

    pair<int, int> bLocalRows = split1D(bSize, processArrayOrder, myRow);
    pair<int, int> bLocalCols = split1D(bSize, processArrayOrder, myCol);

    // Create grids
    int aRows = aLocalRows.second;
    int aCols = aLocalCols.second;
    int bRows = bLocalRows.second;
    int bCols = bLocalCols.second;
    int cRows = aRows;
    int cCols = bCols;

    vector<double> aLocal(aRows * aCols);
    vector<double> bLocal(bRows * bCols);
    vector<double> cLocal(cRows * cCols, 0.0);

    // Index helpers
    auto aAt = [&](int i, int j) -> double& { 
        return aLocal[i * aCols + j]; 
    };
    
    auto bAt = [&](int i, int j) -> double& { 
        return bLocal[i * bCols + j]; 
    };
    
    auto cAt = [&](int i, int j) -> double& { 
        return cLocal[i * cCols + j]; 
    };

    // Fill grids of this process with random values
    for (int i = 0; i < aRows; i++) {
        for (int j = 0; j < aCols; j++) {
            aAt(i, j) = valueAt(i + aLocalRows.first, j + aLocalCols.first, baseSeed);
        }
    }
    // Different base seed to prevent correlation
    // The const is the 'golden ratio'
    uint64_t bSeed = baseSeed ^ 0x9e3779b9;
    for (int i = 0; i < bRows; i++) {
        for (int j = 0; j < bCols; j++) {
            bAt(i, j) = valueAt(i + bLocalRows.first, j + bLocalCols.first, bSeed);
        }
    }

    // Create communicators for broadcasting to row and column 
    MPI_Comm rowComm;
    MPI_Comm colComm;

    MPI_Comm_split(MPI_COMM_WORLD, myRow, myCol, &rowComm);

    MPI_Comm_split(MPI_COMM_WORLD, myCol, myRow, &colComm);

    vector<double> aBcast;
    vector<double> bBcast;


    // Perform multiple runs
    for (int run = 1; run <= runs; run++) {
        // Reset C
        fill(cLocal.begin(), cLocal.end(), 0.0);

        // Synchronize before starting timer for accurate timing
        MPI_Barrier(MPI_COMM_WORLD);
        auto start = chrono::high_resolution_clock::now();

        for (int k = 0; k < processArrayOrder; k++) {
            // Find the size of the broadcast
            pair<int, int> kCols = split1D(aSize, processArrayOrder, k);
            int panelWidth = kCols.second;
            aBcast.resize(aRows * panelWidth);

            if (myCol == k) {
                // copy local A into the buffer
                aBcast = aLocal;
            }
            // Broadcast
            MPI_Bcast(aBcast.data(), aRows * panelWidth, MPI_DOUBLE, k, rowComm);

            bBcast.resize(bCols * panelWidth);

            if (myRow == k) {
                // copy local B into the buffer
                bBcast = bLocal;
            }
            // Broadcast
            MPI_Bcast(bBcast.data(), bRows * panelWidth, MPI_DOUBLE, k, colComm);

            // Perform calculation
            for (int i = 0; i < aRows; i++) {
                for (int p = 0; p < panelWidth; p++) {
                    for (int j = 0; j < bCols; j++) {
                        cAt(i, j) += aBcast[i * panelWidth + p] * bBcast[p * bCols + j];
                    }
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        auto end = chrono::high_resolution_clock::now();

        double localSum = 0.0;
        for (double val : cLocal) {
            localSum += val;
        }

        double globalSum = 0.0;
        MPI_Reduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


        if (rank == 0) {
            cout << "<=========================================================>\n"
                << "RUN: " << run << '\n';
            cout << "Checksum: " << globalSum << "\n";

            double elapsed = chrono::duration<double>(end - start).count();
            cout << "Time: " << elapsed << "s\n";
        }
    }
// P processes
// Each Pij owns Aij, Bij, Cij

// Each process gets its equal split of A and B and C
// runs a for loop
/*
* for k:
    broadcast A column k across rows
    broadcast B row k down columns

    for each process (i,j)
        C(i,j) += A(i,k) * B(k,j)
*/
}



inline uint64_t PM2_DGEMM::Xoshiro256::splitmix64(uint64_t& x) {
    x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

PM2_DGEMM::Xoshiro256::Xoshiro256(uint64_t seed) {
    for (int i = 0; i < 4; i++) {
        s[i] = splitmix64(seed);
    }
}

inline uint64_t PM2_DGEMM::Xoshiro256::rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t PM2_DGEMM::Xoshiro256::next() {
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

uint64_t PM2_DGEMM::Xoshiro256::nextInRange(uint64_t min, uint64_t max) {
    uint64_t range = max - min + 1;
    uint64_t x;
    uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

    do {
        x = next();
    } while (x >= limit);

    return min + (x % range);
}

double PM2_DGEMM::Xoshiro256::nextInRangeDouble(double min, double max) {
    uint64_t x = next();
    const uint64_t mantissaMask = (1ULL << 53) - 1;

    double zto = (x >> 11) * (1.0 / (1ULL << 53));
    return min + (max - min) * zto;
}
