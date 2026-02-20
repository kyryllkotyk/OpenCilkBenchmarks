#ifndef PM2_DGEMM_
#define PM2_DGEMM_

#include <vector>
#include <iostream>
using namespace std;

class PM2_DGEMM
{
public:
    /*
    @brief Performs tiled matrix multiplication C = A × B.

    @param runs How many runs to perform
    @param A First matrix for multiplication
    @param B Second matrix for multiplication
    @param blockRowSize Row size of each block (submatrix of C)
    @param blockColSize Column size of each block (submatrix of C)
    @param baseSeed RNG seed to start from (incremented by rank)

    @return Result A x B matrix multiplication. Size Arow x Bcol.
    Returns empty matrix on failure / misinput
    */
    vector<vector<double>> runBenchmark(const short runs,
        /*Distributed Info*/int rank,
        const vector<vector<double>>& A, const vector<vector<double>>& B,
        const short blockRowSize, const short blockColSize, int baseSeed);

private:
    struct Xoshiro256 {
    private:
        uint64_t s[4];

        static inline uint64_t splitmix64(uint64_t& x);
        static inline uint64_t rotl(uint64_t x, int k);

    public:
        explicit Xoshiro256(uint64_t seed);
        uint64_t next();
        uint64_t nextInRange(uint64_t min, uint64_t max);
        double nextInRangeDouble(double min, double max);
    };
};

#endif
