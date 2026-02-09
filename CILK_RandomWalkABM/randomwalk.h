#ifndef RANDOM_WALK_
#define RANDOM_WALK_

#include <vector>
#include <algorithm> // For swap
#include <string>
#include <limits>
#include <iostream>
using namespace std;

class RandomWalk
{
private:
  struct Agent {
    Agent(int id, int xI, int yI, int wealthI, int metabolismI, int visionI)
      : ID(id), x(xI), y(yI), wealth(wealthI), metabolism(metabolismI)
      , vision(visionI) {};

    Agent() = default;

    int ID;
    int x, y;
    int wealth;
    int metabolism;
    int vision;
    /*
    // Components (optional)
    int age = 0;
    int deathAge = 100;
    // True = Male, False = Female
    bool maleSex = true;
    */
  };

  struct GridCell {
    GridCell() = default;

    GridCell(int capacity) : sugarCapacity(capacity), currentSugar(capacity) {};

    int sugarCapacity;
    int currentSugar;
  };

  struct Xoshiro256 {
  private:
    uint64_t s[4];

    static inline uint64_t splitmix64(uint64_t& x);
    static inline uint64_t rotl(uint64_t x, int k);

  public:
    explicit Xoshiro256(uint64_t seed);    
    uint64_t next();
    uint64_t nextInRange(uint64_t min, uint64_t max);
  };


public:
  RandomWalk() = default;
  //TODO:: add a flag for interactions to enable or disable
  /*
  * @param grid Reference to the grid that agents walk on. 
    Cell value is the agent count. Will be modified by the function
  * @param agents Collection of agents to walk on the grid
  * @param timeSteps How many timesteps to simulate
  */
  void runBenchmarkPrototype(vector<vector<int>>& grid, 
    vector<RandomWalk::Agent>& agents, const short timeSteps);

  /*
  struct ModuleParameters {
    unsigned int minDeathAge = UINT_MAX;
    unsigned int maxDeathAge = UINT_MAX;
    unsigned int fertilityMin = UINT_MAX;
    unsigned int fertilityMax = UINT_MAX;
    unsigned int reproductionThreshold = UINT_MAX;
    double inheritanceFraction = numeric_limits<double>::quiet_NaN();
  };
  */

  //TODO:: write out description for inputs
  // for models to use -> "ar"/"ra" = reproduction + aging, "a" = aging, "r" = reproduction
  void runBenchmark(const unsigned int gridX, const unsigned int gridY, 
    const unsigned int agentCount,
    const unsigned int growthRate, const unsigned int timeSteps,
    const unsigned int sugarCapacityMin, const unsigned int sugarCapacityMax,
    const unsigned int metabolismMin, const unsigned int metabolismMax,
    const unsigned int visionMin, const unsigned int visionMax,
    unsigned int initialWealth = UINT_MAX
    /*
    const string& modulesToUse = "",
    const ModuleParameters& moduleParameters = {}
    */
   );

  

  // SUGARSCAPE
  // Agent:
  //  ID
  //  Location (randomized at start)
  //  Wealth (how much sugar the agent has)
  //  Metabolism (per timestep survival cost, so wealth -= metabolism each timestep)
  //   
  // Constants:
  //  Growth Rate (How much sugar grows in a cell after each timestep)
  //  Total time steps
  //  Initial wealth per-agent (default to between metabolism * 5 - * 10)
  //    Randomization Seed: 100
  //  Sugar Capacity Min (minimum sugar capacity a grid can have)
  //  Sugar Capacity Max (maximum sugar capacity a grid can have)
  //  Metabolism Min (minimum metabolism value)
  //  Metabolism Max (maximum metabolism value)
  //  Vision Min (minimum vision value)
  //  Vision Max (maximum vision value)
  // 
  // Grid Cells:
  //  Sugar Count (how much sugar is currently present)
  //    Starts with sugar capacity at t = 0
  //  Sugar Capacity (max amount of sugar in that cell)
  //
  // Notes:
  //  Sugar is regrown in the beginning of a step
  //  Allow staying in the same cell, force move only if better cell exists
  //    Agent's current cell is a valid candidate
  //  Phase 1 - agents choose target cell WITHOUT moving 
  //    At equal sugar count, move to the closest cell
  //      At equal distance (measured via steps), choose the lowest (x,y), with x preference
  //  Phase 2 - conflict resolution, winner moves loser stays
  //    Lowest ID value wins
  //  Agent takes all the sugar out of the cell
  //  Fixed seed pre-generated sugar capacity via Xoshiro256
  //    Seed: 1
  //  Fixed seed pre-generated metabolism via Xoshiro256
  //    Seed: 2
  //  Fixed seed pre-generated vision via Xoshiro256
  //    Seed: 3
  //  Fixed seed pre-generated coordinates via Xoshiro256
  //    Seed: 4
  // 
  // Canonical Modules:
  //  Aging -
  //    Agents:
  //      Age (timesteps it has lived so far)
  //    Constants:
  //      Age Min (minimum death age that can be assigned to agent)
  //      Age Max (maximum death age that can be assigned to agent)
  //    Notes:
  //      Fixed pre-generated age via Xoshiro256
  //          Seed: 5
  //      When the number of timesteps it has been alive reaches death age,
  //      the agent dies, regardless of wealth
  //
  //  Reproduction -
  //    Agent:
  //      Sex (Male/Female)
  //    Constants:
  //      Fertility Age Min (minimum age at which reproduction can happen)
  //      Fertility Age Max (maximum age at which reproduction can happen)
  //      Reproduction Threshold (minimum wealth to reproduce)
  //      Inheritance fraction (what percentage of each parent's wealth the child gets)
  //    Notes:
  //      To reproduce, the agent must be
  //        Within the fertility age range
  //        At least one adjacent agent (4 way cardinal neighbors)
  //          Opposite sex
  //          Fertile
  //          Also above wealth threshold
  //        At least one empty adjacent cell (to the current cell, not other parent)
  //      When agents reproduce
  //        A new agent is created in an adjacent cell
  //          Lowest (x,y) adjacent empty cell
  //        Both parents transfer wealth to the child
  //           Each parent gives inheritance fraction of their respective wealth
  //        Child gets
  //          Age = 0
  //          Wealth = Inherited Amount
  //          Randomized vision via same Xoshiro256 instance as pre-generation
  //              Seed: 3 (Re-use the vision generator, don't create a new one)
  //          Randomized sex via seeded Xoshiro256
  //            Seed: 6
  //      Reproduction must only happen once per timestep per agent, in increasing ID order



  // Potentially change phase 2 to make the loser pick the second best spot
};

#endif
