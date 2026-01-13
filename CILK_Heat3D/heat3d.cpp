#include "heat3d.h"

vector<pair<uint64_t, uint64_t>> Heat3D::runBenchmark(vector<float>& grid,
  const short gridX, const short gridY, const short gridZ,
  const short timesteps, const float alpha, const float beta,
  const short decompX, const short decompY, const short decompZ, 
  const short totalRuns, const short warmupRuns) {
  
  // Ensure that the grid size matches the given size parameters
  if (grid.size() != (gridX + 2) * (gridY + 2) * (gridZ + 2)) {
    exit(0);
  }
  
  // Write down the padded sizing to not recalculate over and over
  int Y = gridY + 2;
  int Z = gridZ + 2;
  // Number of elements to skip when moving +1 in x
  int strideX = Y * Z;

  // index formula = x * (Y * Z) + y * Z + z
  // where x,y,z = for loop values
  // Y, Z = gridY + 2, gridZ + 2
#define IDX(x,y,z) ((x) * strideX + (y) * Z + (z))

  // Storage for halo & execution times
  vector<pair<uint64_t, uint64_t>> timesPerRunUs(totalRuns);
  // Accumulated halo & execution times PER RUN
  uint64_t accumulatedHaloTime = 0, accumulatedExecutionTime = 0;

  // Grid for t + 1
  vector<float> newGrid = grid;

  // Run the entire benchmark multiple times
  for (int runs = 0; runs < totalRuns + warmupRuns; runs++) {
    // Simulate specified timesteps
    for (int i = 0; i < timesteps; i++) {

      // Start measuring halo time here
      auto haloTimeStart = std::chrono::steady_clock::now();

      // Halos - Set the boundaries to be equal to the nearest valid cell 
      // on the respective axis 
      // X halos
      for (int y = 1; y <= gridY; y++) {
        for (int z = 1; z <= gridZ; z++) {
          grid[IDX(0, y, z)] = grid[IDX(1, y, z)];
          grid[IDX(gridX + 1, y, z)] = grid[IDX(gridX, y, z)];
        }
      }
      // Y halos
      for (int x = 0; x <= gridX + 1; x++) {
        for (int z = 1; z <= gridZ; z++) {
          grid[IDX(x, 0, z)] = grid[IDX(x, 1, z)];
          grid[IDX(x, gridY + 1, z)] = grid[IDX(x, gridY, z)];
        }
      }
      // Z halos
      for (int x = 0; x <= gridX + 1; x++) {
        for (int y = 0; y <= gridY + 1; y++) {
          grid[IDX(x, y, 0)] = grid[IDX(x, y, 1)];
          grid[IDX(x, y, gridZ + 1)] = grid[IDX(x, y, gridZ)];
        }
      }

      // Stop measuring halo time here
      auto haloTimeEnd = std::chrono::steady_clock::now();
      
      // Start measuring execution time here
      auto executionTimeStart = std::chrono::steady_clock::now();
      // Loop over the blocks
      cilk_for (int blockId = 0; blockId < decompX * decompY * decompZ; blockId++) {
        // Map 1D index to 3D block coordinates (bz varies fastest)
        int bz = blockId % decompZ;
        int by = (blockId / decompZ) % decompY;
        int bx = blockId / (decompZ * decompY);

        // Find the range of this block in the grid, including the offset
        int x0 = (bx * gridX) / decompX + 1;
        int x1 = (bx + 1) * gridX / decompX + 1;

        int y0 = (by * gridY) / decompY + 1;
        int y1 = ((by + 1) * gridY) / decompY + 1;

        int z0 = (bz * gridZ) / decompZ + 1;
        int z1 = ((bz + 1) * gridZ) / decompZ + 1;

        // Loop over all the cells in the block
        for (int x = x0; x < x1; x++) {
          for (int y = y0; y < y1; y++) {
            for (int z = z0; z < z1; z++) {
              // index of the current cell
              int c = IDX(x, y, z);
              // Apply the update formula
              newGrid[c] = grid[c] * alpha + beta *
                (grid[IDX(x - 1, y, z)] + grid[IDX(x + 1, y, z)] +
                  grid[IDX(x, y - 1, z)] + grid[IDX(x, y + 1, z)] +
                  grid[IDX(x, y, z - 1)] + grid[IDX(x, y, z + 1)]);
            }
          }
        }
      }
      // Stop measuring execution time here
      auto executionTimeEnd = std::chrono::steady_clock::now();
      // Since T is becoming Tprev + 1, Tprev + 1 values are new T values 
      swap(grid, newGrid);
      // Accumulate halo time from this loop
      accumulatedHaloTime += std::chrono::duration_cast
        <std::chrono::microseconds>(haloTimeEnd -
          haloTimeStart).count();
      // Accumulate execution time from this loop 
      accumulatedExecutionTime += std::chrono::duration_cast
        <std::chrono::microseconds>(executionTimeEnd -
          executionTimeStart).count();
    }

    // Record times only if it's no longer a warmup run
    if (runs >= warmupRuns) {
      timesPerRunUs[runs - warmupRuns] = { accumulatedHaloTime, accumulatedExecutionTime };
    }
    accumulatedHaloTime = 0;
    accumulatedExecutionTime = 0;
  }

  return timesPerRunUs;
}
