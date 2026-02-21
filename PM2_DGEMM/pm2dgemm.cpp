#include "pm2dgemm.h"

inline pair<int,int> split1D(int totalDimensionSize, int processesInDimension, int index) {
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
    vector<vector<double>> aLocalGrid(aLocalRows.second, vector<double>(aLocalCols.second));
    vector<vector<double>> bLocalGrid(bLocalRows.second, vector<double>(bLocalCols.second));
    vector<vector<double>> cLocalGrid(aLocalRows.second, vector<double>(bLocalCols.second));

    // Fill grid of this process with random values
    for (int i = 0; i < aLocalRows.second; i++) {
        for (int j = 0; j < aLocalCols.second; j++) {
            aLocalGrid[i][j] = valueAt(i + aLocalRows.first, j + aLocalCols.first, baseSeed);
        }
    }
    // Different base seed to prevent correlation
    // The const is the 'golden ratio'
    uint64_t bSeed = baseSeed ^ 0x9e3779b9;
    for (int i = 0; i < bLocalRows.second; i++) {
        for (int j = 0; j < bLocalCols.second; j++) {
            bLocalGrid[i][j] = valueAt(i + bLocalRows.first,
                j + bLocalCols.first, bSeed);
        }
    }

    for (int k = 0; k < processArrayOrder; k++) {
        if (k == myCol) 
            //broadcast A for processes in the same row
            // Buffer, count, MPI datatype, root, 
            MPI_Bcast(aLocalGrid[0], aLocalRows, rank, MPI_DOUBLE, )
        }
        else {
            //receive A from process[myRow][k]
        }

        if (k == myRow){
            //broadcast B for processes in the same column
        }
        else {
            //receive B from process[k][myCol]
        }

        //Perform calculation
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
    return vector<vector<double>>();
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
