#include "heat3d.h"

// forward declaration
static std::string pad5(int x); 

vector<pair<uint64_t, uint64_t>> Heat3D::runBenchmark(vector<double>& grid,
  const short gridX, const short gridY, const short gridZ,
  const short timesteps, const float alpha, const float beta,
  const short decompX, const short decompY, const short decompZ,
  const short totalRuns, const short warmupRuns,
  const short screenshotEvery, const char borderType) {
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
  vector<double> newGrid = grid;

  // Run the entire benchmark multiple times
  for (int runs = 0; runs < totalRuns + warmupRuns; runs++) {

    std::vector<int> writtenTimesteps;

    // Simulate specified timesteps
    for (int i = 0; i < timesteps; i++) {
      // Save the grid at the beginning of the step 
      // Grid, timestep, locality (const 1)
      if (runs == 0 && screenshotEvery != 0 && i % screenshotEvery == 0) {
        writeToFile(grid, i, 0, gridX, gridY, gridZ, "./out_vtk/");
        writtenTimesteps.push_back(i);
      }

      switch (borderType) {
      case('d'): /* Dirichlet */
        // Start measuring halo time here
        auto haloTimeStart = std::chrono::steady_clock::now();

        // Halos - Set the boundaries to be equal to the nearest valid cell 
        // on the respective axis 
        // X halos
        for (int y = 1; y <= gridY; y++) {
          for (int z = 1; z <= gridZ; z++) {
            grid[IDX(0, y, z)] = 0;
            grid[IDX(gridX + 1, y, z)] = 0;
          }
        }
        // Y halos
        for (int x = 0; x <= gridX + 1; x++) {
          for (int z = 1; z <= gridZ; z++) {
            grid[IDX(x, 0, z)] = 0;
            grid[IDX(x, gridY + 1, z)] = 0;
          }
        }
        // Z halos
        for (int x = 0; x <= gridX + 1; x++) {
          for (int y = 0; y <= gridY + 1; y++) {
            grid[IDX(x, y, 0)] = 0;
            grid[IDX(x, y, gridZ + 1)] = 0;
          }
        }

      case('i'): /* Insulated */
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

      case('h'): /* Heater Window */
        break;
      }

      // Stop measuring halo time here
      auto haloTimeEnd = std::chrono::steady_clock::now();

      // Start measuring execution time here
      auto executionTimeStart = std::chrono::steady_clock::now();
      // Loop over the blocks
      cilk_for(int blockId = 0; blockId < decompX * decompY * decompZ; blockId++) {
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

    // Final state
    if (runs == 0 && screenshotEvery > 0) {
      writeToFile(grid, timesteps, 0, gridX, gridY, gridZ, "./out_vtk/");
      writtenTimesteps.push_back(timesteps);
      writePvd(writtenTimesteps, 1, "./out_vtk/");
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





void Heat3D::writeToFile(vector<double>& grid, int timestep, int locality,
  int nx, int ny, int nz, string directory)
{
  if (grid.size() != (nx + 2) * (ny + 2) * (nz + 2)) {
    return;
  }

  // Create directory (if needed)
  int rc = mkdir(directory.c_str(), 0755);
  if (rc != 0 && errno != EEXIST) {
    throw std::runtime_error(
      "mkdir failed for '" + directory + "': " + std::string(strerror(errno)));
  }

  // Generate the file name
  // {directory}/heat3d_t{timestep:05d}_loc{locality:05d}.vti
  std::string fileName =
    directory + "/heat3d_t" + pad5(timestep) + "_loc" + pad5(locality) + ".vti";

  // Calculate extents (for cilk, whole grid (0 to (n - 1), so skip)
  // Write the XML header structure
  std::ofstream out(fileName.c_str());
  if (!out.is_open()) {
    throw std::runtime_error("Cannot open file: " + fileName);
  }

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <ImageData WholeExtent=\""
    << 0 << " " << nx - 1 << " "
    << 0 << " " << ny - 1 << " "
    << 0 << " " << nz - 1 << "\""
    << " Origin=\"0.0 0.0 0.0\""
    << " Spacing=\"1.0 1.0 1.0\">\n";

  out << "    <Piece Extent=\""
    << 0 << " " << nx - 1 << " "
    << 0 << " " << ny - 1 << " "
    << 0 << " " << nz - 1 << "\">\n";

  out << "      <PointData Scalars=\"temperature\">\n";
  out << "        <DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\" NumberOfComponents=\"1\">\n";
  out << std::scientific << std::setprecision(15);

  // Write the temperature values

  // haloed dimensions
  int Y = ny + 2;
  int Z = nz + 2;
  int strideX = Y * Z;

  // Write interior values only, in VTK order (k, j, i)
  for (int k = 0; k < nz; k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {

        // shift into haloed storage
        int si = i + 1;
        int sj = j + 1;
        int sk = k + 1;

        // 1D index into haloed vector (row-major)
        size_t idx = static_cast<size_t>(si) * static_cast<size_t>(strideX)
          + static_cast<size_t>(sj) * static_cast<size_t>(Z)
          + static_cast<size_t>(sk);

        out << grid[idx];

        // space between values (avoid trailing space at very end)
        if (!(i == nx - 1 && j == ny - 1 && k == nz - 1)) {
          out << " ";
        }
      }
    }
  }
  out << "\n";

  // Write XML footer and close file
  out << "        </DataArray>\n";
  out << "      </PointData>\n";
  out << "      <CellData>\n";
  out << "      </CellData>\n";
  out << "    </Piece>\n";
  out << "  </ImageData>\n";
  out << "</VTKFile>\n";
  out.close();

}

static std::string pad5(int x) {
  std::ostringstream ss;
  ss << std::setw(5) << std::setfill('0') << x;
  return ss.str();
}

void Heat3D::writePvd(const std::vector<int>& writtenTimesteps,
  int numLocalities, const std::string& directory) {
  if (writtenTimesteps.empty()) return;
  if (numLocalities <= 0) return;

  std::string pvdName = directory + "/heat3d.pvd";
  std::ofstream out(pvdName.c_str());
  if (!out.is_open()) {
    throw std::runtime_error("Cannot open file: " + pvdName);
  }

  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <Collection>\n";

  for (size_t ti = 0; ti < writtenTimesteps.size(); ++ti) {
    int t = writtenTimesteps[ti];

    for (int loc = 0; loc < numLocalities; ++loc) {
      // IMPORTANT: relative path only (just filename)
      std::string vtiName = "heat3d_t" + pad5(t) + "_loc" + pad5(loc) + ".vti";

      out << "    <DataSet timestep=\"" << t
        << "\" group=\"\" part=\"" << loc
        << "\" file=\"" << vtiName << "\"/>\n";
    }
  }

  out << "  </Collection>\n";
  out << "</VTKFile>\n";
  out.close();
}
