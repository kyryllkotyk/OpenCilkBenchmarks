#include "dgemm.h"

vector<vector<double>> DGEMM::runBenchmark(
  const vector<vector<double>>& A, const vector<vector<double>>& B,
  const short blockRowSize, const short blockColSize) {
  
  // Output vector
  vector<vector<double>> C;

  // Small error handling to ensure valid input (in scope)
  if (A.size() == 0 || B.size() == 0 || A[0].size() == 0 || B[0].size() == 0) {
    return C;
  }
  if (blockRowSize <= 0 || blockColSize <= 0) {
    return C;
  }

  // Ensure that the matrices CAN be multiplied 
  // AxB x BxC = AxC
  if (A[0].size() != B.size()) {
    return C;
  }

  // Number of rows in A
  int ra = A.size();
  // Number of columns in A & rows in B
  int carb = A[0].size();
  // Number of columns in B and C
  int cbc = B[0].size();

  // AxC result vector to update
  C = vector<vector<double>>(ra, vector<double>(cbc, 0));

  // Find how many tiles/blocks we need to cover the output matrix C
  // If ra or cbc is a multiple of the block size, it is rounded up
  // so there's a partial tile for leftovers
  int tileRows = (ra + blockRowSize - 1) / blockRowSize;
  int tileCols = (cbc + blockColSize - 1) / blockColSize;
  int totalBlocks = tileRows * tileCols;

  // Declared here as to not redeclare it inside the for loop each time
  int iEnd, jEnd; 


  // Iterate over every output tile with a single id
  // Makes it easy to paralellize with cilk_for as each blockId corresponds to
  // a subspace of C to calculate
  // As long as blockRow and blockCol are unique, the region of C is unique,
  // so no race conditions.
  // TODO:: replace with Cilk for, add timers and other data gathering
  cilk_for (int blockId = 0; blockId < totalBlocks; blockId++) {
    // Convert blockId into 2D tile coordinates
    // blockRow is the which tile row it's in
    // blockCol is the tile column it's in
    // blockId = blockRow * tileCols + blockCol
    int blockRow = blockId / tileCols;
    int blockCol = blockId % tileCols;

    // Starting (top left) indices in C for this tile
    // Each tile covers i/j to i/j + blockRow/ColSize - 1
    int i = blockRow * blockRowSize;
    int j = blockCol * blockColSize;

    // Compute the exclusive end indices for this tile, clipped to the 
    // matrix size so edge tiles don't go out of bounds
    // iEnd is at most ra (number of rows in A/C)
    // jEnd is at most cbc (number of cols in B/C)
    iEnd = std::min(i + (int)blockRowSize, ra);
    jEnd = std::min(j + (int)blockColSize, cbc);

    // This order (ii -> k -> jj) is chosen because 
    // each value from A is loaded once and reused to update 
    // multiple columns of C instead of reloading it for every column
    // For this tile of C:
    // for each row ii in the tile
    //  for each shared dimension k
    //    for each col jj in the tile
    //      C[ii][jj] += A[ii][k] * B[k][jj]
    for (int ii = i; ii < iEnd; ii++) {
      for (int k = 0; k < carb; k++) {
        double aVal = A[ii][k];
        for (int jj = j; jj < jEnd; jj++) {
          C[ii][jj] += aVal * B[k][jj];
        }
      }
    }
  }

  return C;
}
