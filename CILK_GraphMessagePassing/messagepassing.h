#ifndef MESSAGE_PASSING_
#define MESSAGE_PASSING_

#include <vector>

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

using namespace std;

class MessagePassing
{
public:
	MessagePassing() = default;
	//TODO:: modify vertices array to adjust to more complicated messages, if needed / Ahmed decides.
	void runSimulation(vector<vector<int>>& adjList, vector<int>& vertices);
};

#endif
