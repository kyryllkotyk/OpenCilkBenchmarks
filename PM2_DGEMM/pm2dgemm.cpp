#include "pm2dgemm.h"

vector<vector<double>> PM2_DGEMM::runBenchmark(const short runs,
    int aRowFullSize, int aColFullSize, int bRowFullSize, int bColFullSize,
    int blockRowSize, int blockColSize,
    int baseSeed) {

    if (aColFullSize != bRowFullSize) {
        cerr << "Cannot multiply, expected sizes: XxY, YxZ";
        return {};
    }

    // Figure out the rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Figure out how many ranks there are
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Seed the randomizer for this process
    // For uniqueness, add rank to the base seed
    // Otherwise, same seed -> same values
    Xoshiro256 rng(baseSeed + rank);

    // Turn the process list into a 2D array
    int dims[2] = { 0, 0 };
    MPI_Dims_create(size, 2, dims);
    int processesRow = dims[0], processesCol = dims[1];

    // Figure out grid dimensions
    int myRow = rank / processesCol;
    int myCol = rank % processesCol;

    // Create the A,B,C matrices for this process
    // blockSize x blockSize
    vector<vector<double>> A, B, C;

    // Fill grid of this process
    
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
