#ifndef RANDOM_WALK_
#define RANDOM_WALK_

#include <vector>
#include <algorithm> // For swap

class RandomWalk
{
private:
  struct Agent {
    //TODO::implement
    int setCoordinates(int xS, int yS);
    int x = 0, y = 0; // Coordinates
    //TODO:: Add optional state here for interaction
    //If NOT adding optional state, then remove Agent and
    // refactor the vector in runBenchmark to be a pair of x,y coordinates
  };

  struct Xoshiro256 {
    uint64_t s[4];

    static inline uint64_t splitmix64(uint64_t& x);

    explicit Xoshiro256(uint64_t seed);

    static inline uint64_t rotl(uint64_t x, int k);

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
  void runBenchmark(std::vector<std::vector<int>>& grid, 
    std::vector<RandomWalk::Agent>& agents, const short timeSteps);
};

#endif
