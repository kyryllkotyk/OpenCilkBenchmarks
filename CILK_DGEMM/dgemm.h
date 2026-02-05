#ifndef DGEMM_
#define DGEMM_

#include <vector>

#include "cilk/cilk.h"
using namespace std;

class DGEMM
{
public:
 
  DGEMM() = default;
  /*
    Performs tiled matrix multiplication C = A × B.

    @param A First matrix for multiplication
    @param B Second matrix for multiplication
    @param blockRowSize Row size of each block (submatrix of C)
    @param blockColSize Column size of each block (submatrix of C)
    @return Result A x B matrix multiplication. Size Arow x Bcol.
    Returns empty matrix on failure / misinput
  */
  vector<vector<double>> runBenchmark(
    const vector<vector<double>>& A, const vector<vector<double>>& B,
    const short blockRowSize, const short blockColSize);

};

#endif
