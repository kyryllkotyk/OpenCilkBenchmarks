#ifndef HEAT_3D_
#define HEAT_3D_

#include <vector>
#include <algorithm> // For swap

#include <chrono> // For timing

#include <cilk/cilk.h> // CILK

#include <sys/stat.h>  // mkdir
#include <sys/types.h>
#include <cerrno>      // errno

#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstring>
#include <stdexcept>

using namespace std;
class Heat3D
{
public:
  Heat3D() = default;

  /*
    @brief Splits the given grid into blocks and calculates the heat value for every cell.
    Modifies the given grid. 

    @param grid 1D array with the size of (gridX + 2) * (gridY + 2) * (gridZ + 2).
    To convert from a 3D array, perform row flattening with z as the fast varying
    index, INCLUDING the halo padding (0th and dimension size # elements)
    @param gridX, gridY, gridZ Dimension size of the grid (Without padding)
    @param timesteps How many steps (updates) to perform
    @param alpha Constant value for the update formula for this cell's value importance
    @param beta Constant value for the update formula for neighbor cell's value importance
    @param decompX, decompY, decompZ Dimension size of the block to split grid into 
    (Without padding)
    @param totalRuns How many times to perform the same benchmark
    @param warmupRuns How many times to perform the benchmark BEFORE recording
    @param borderType 'd' - Dirichlet, 'i' - Insulated, 'h' - Heater Window

    @return Halo time & execution time for each run
  */
  vector<pair<uint64_t, uint64_t>> runBenchmark(vector<double>& grid, 
    const short gridX, const short gridY, const short gridZ, 
    const short timesteps, const float alpha, const float beta,
    const short decompX, const short decompY, const short decompZ,
    const short totalRuns, const short warmupRuns, const short screenshotEvery,
    const char borderType);

  // Writes to a VTI file based on specified file format:
  // heat3d_t{timestep:05d}_loc{locality:05d}.vti
  void writeToFile(vector<double>& grid, int timestep, int locality,
    int nx, int ny, int nz, string directory);

  void writePvd(const std::vector<int>& writtenTimesteps,
    int numLocalities, const std::string& directory);
};

#endif
