#include "randomwalk.h"

void RandomWalk::runBenchmarkPrototype(std::vector<std::vector<int>>& grid,
  std::vector<Agent>& agents, const short timeSteps) {
  
  //TODO:: Add guarding for misinput
  
  // Deterministic seeded randomizer
  Xoshiro256 randomizer(1000);

  int gridX = grid.size();
  int gridY = grid[0].size();

  std::vector<Agent> nextAgents = agents;
  /*
  State at each time step:
    Set of agents
      Unique ID (Vector index in this case)
      Position (X,Y)
      (OPTIONAL) Internal State
    Spacial Domain
      2D grid sized Lx x Ly

  FOR EACH TIMESTEP:
  */
  for (int step = 0; step < timeSteps; step++) {
    /*
    Phase 1: Movement Direction
      Action:
        For every agent
          Read current position
          Get a random number
            ** USING MOORE MOVEMENT DIRECTION (8 way) **
            Either
              1 - 8 signifying direction
              +- 1 for x, +- 1 for y
              1 - 3 for each, x and y, where 1 = -1, 2 = 0, 3 = 1
              2nd option is likely to be used
          Find the agent's next position based on the random number
          Write that position down in the t + 1 agent list
    */
    for (int i = 0; i < agents.size(); i++) {
      // 1 = Move to the left/down
      // 2 = Don't move
      // 3 = Move to the right/up
      int rand = randomizer.nextInRange(1, 3);
      int newX = agents[i].x + rand - 2;
      // Clamp X
      if (newX >= gridX) {
        newX = gridX - 1;
      }
      if (newX < 0) {
        newX = 0;
      }
      nextAgents[i].x = newX;

      // Repeat for y
      rand = randomizer.nextInRange(1, 3);
      int newY = agents[i].y + rand - 2;
      // Clamp Y
      if (newY >= gridY) {
        newY = gridY - 1;
      }
      if (newY < 0) {
        newY = 0;
      }
      nextAgents[i].y = newY;

    }

    /*
    Phase 2: Movement Resolution
        Action:
          Execute the move based to the proposed position from phase 1
            Swap the current agent list with agent list of t + 1
          Rebuild grid from the new agent list
    */
    std::swap(agents, nextAgents);

    // Clear grid
    for (auto& col : grid) {
      std::fill(col.begin(), col.end(), 0);
    }

    // Fill the grid with updated agent locations
    for (Agent& agent : agents) {
      grid[agent.x][agent.y]++;
    }

    /*
    Phase 3 (OPTIONAL): Interactions
        Action:
          For each agent
            Identify their neighbor (agents in same or adjacent cells)
            Do the interaction
            Update agent's state (NOT THE NEIGHBOR'S)
    */
    //TODO:: Add interactions if needed
  }

}

inline uint64_t RandomWalk::Xoshiro256::splitmix64(uint64_t& x) {
  x += 0x9e3779b97f4a7c15ULL;
  uint64_t z = x;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

RandomWalk::Xoshiro256::Xoshiro256(uint64_t seed) {
  for (int i = 0; i < 4; i++) {
    s[i] = splitmix64(seed);
  }
}

inline uint64_t RandomWalk::Xoshiro256::rotl(uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

uint64_t RandomWalk::Xoshiro256::next() {
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

uint64_t RandomWalk::Xoshiro256::nextInRange(uint64_t min, uint64_t max) {
  uint64_t range = max - min + 1;
  uint64_t x;
  uint64_t limit = UINT64_MAX - (UINT64_MAX % range);

  do {
    x = next();
  } while (x >= limit);

  return min + (x % range);
}


void RandomWalk::runBenchmark(
  const unsigned int height, const unsigned int width, const unsigned int agentCount,
  const unsigned int growthRate, const unsigned int timeSteps,
  const unsigned int sugarCapacityMin, const unsigned int sugarCapacityMax,
  const unsigned int metabolismMin, const unsigned int metabolismMax,
  const unsigned int visionMin, const unsigned int visionMax, 
  unsigned int initialWealth) {
  
  // Error checking block
  // Empty grid
  if (height == 0 || width == 0) {
    cerr << "Height or width set to 0\n";
    return;
  }

  // Too many agents
  if (height * width < agentCount) {
    cerr << "Too many agents to spawn";
    return;
  }
  
  // Invalid sugar capacity 
  if (sugarCapacityMax < sugarCapacityMin) {
    cerr << "Invalid sugar capacity (Max must be greater or equal to minimum)";
    return;
  }

  // Invalid metabolism
  if (metabolismMax < metabolismMin) {
    cerr << "Invalid metabolism (Max must be greater or equal to minimum)";
    return;
  }

  // Invalid vision
  if (visionMax < visionMin) {
    cerr << "Invalid vision (Max must be greater or equal to minimum)";
    return;
  }

  // Degenerate simulation errors
  if (sugarCapacityMax == 0) {
    cerr << "Degenerate Simulation: No resources possible";
    return;
  }

  if (initialWealth == 0) {
    cerr << "Degenerate Simulation: Guaranteed instant agent death";
    return;
  }

  if (agentCount == 0) {
    cerr << "Degenerate Simulation: No operations possible";
    return;
  }

  /* PRE SIMULATION GRID HANDLING */
  // Initialize sugar capacity randomizer
  Xoshiro256 sugarCapacityRandom(1);
  // Create grid
  vector<vector<GridCell>> grid(height, vector<GridCell>(width));
  // Fill in the grid cell's values using the randomizer
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      grid[i][j].sugarCapacity = grid[i][j].currentSugar = 
        sugarCapacityRandom.nextInRange(sugarCapacityMin, sugarCapacityMax);
    }
  }
  
  /* PRE SIMULATION AGENT HANDLING */
  // Create agent vector
  vector<Agent> agents(agentCount);
  // Initialize metabolism randomizer
  Xoshiro256 metabolismRandom(2);
  // Initialize vision randomizer
  Xoshiro256 visionRandom(3);
  // Initialize wealth randomizer (may be unused)
  Xoshiro256 wealthRandom(100);

  for (int i = 0; i < agentCount; i++) {
    // Give each agent their ID while we're at it
    agents[i].ID = i;
    int meta;
    agents[i].metabolism = meta = metabolismRandom.nextInRange(metabolismMin, metabolismMax);
    agents[i].vision = visionRandom.nextInRange(visionMin, visionMax);
    agents[i].wealth = (initialWealth == UINT_MAX) ? 
      wealthRandom.nextInRange(meta * 5, meta * 10) : initialWealth;
  }

  // To handle coordinates
  // Initialize coordinate randomizer
  Xoshiro256 coordinateRandom(4);
  // Create a temporary vector of all coordinates
  {
    vector<pair<int, int>> allCells;
    allCells.reserve(height * width);

    // Fill the coordinate vector with all possible coordinates
    for (int i = 0; i < height; i++) {
      for (int j = 0; j < width; j++) {
        allCells.push_back({i, j});
      }
    }

    // Shuffle the vector with the coordinate randomizer
    std::shuffle(allCells.begin(), allCells.end(), coordinateRandom);
    
    // Assign the first agentCount # of coordinates to the agents in ascending order
    for (int i = 0; i < agentCount; i++) {
      agents[i].x = allCells[i].first;
      agents[i].y = allCells[i].second;
    }
  }


}
