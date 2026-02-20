#include "pm2dgemm.h"

vector<vector<double>> PM2_DGEMM::runBenchmark(const short runs, 
    int machines, int ranks, const vector<vector<double>>& A,
    const vector<vector<double>>& B, const short blockRowSize, 
    const short blockColSize)
{

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
