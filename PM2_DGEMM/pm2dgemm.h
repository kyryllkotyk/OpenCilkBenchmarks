#ifndef PM2_DGEMM_
#define PM2_DGEMM_

#include <vector>

using namespace std;

class PM2_DGEMM
{
  /*
  @brief Performs tiled matrix multiplication C = A × B.

  @param runs How many runs to perform
  @param A First matrix for multiplication
  @param B Second matrix for multiplication
  @param blockRowSize Row size of each block (submatrix of C)
  @param blockColSize Column size of each block (submatrix of C)
  @return Result A x B matrix multiplication. Size Arow x Bcol.
  Returns empty matrix on failure / misinput
  */
   vector<vector<double>> runBenchmark(const short runs,
       /*Distributed Info*/int machines, int ranks,
       const vector<vector<double>>& A, const vector<vector<double>>& B,
       const short blockRowSize, const short blockColSize);
};

#endif
