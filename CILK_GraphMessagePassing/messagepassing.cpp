#include "messagepassing.h"

void MessagePassing::runSimulation(vector<vector<int>>& adjList,
	vector<int>& vertices, const short timesteps)
{

	int nworkers = __cilkrts_get_nworkers();
	vector<vector<int>> localAccumulator(nworkers, 
		vector<int>(vertices.size(), 0));

	for (int step = 0; step < timesteps; step++) {
		cilk_for (int i = 0; i < adjList.size(); i++) {
			int wid = __cilkrts_get_worker_number();
			for (int j = 0; j < adjList[i].size(); j++) {
				localAccumulator[wid][adjList[i][j]] += 1; //TODO:: change this if a more complicated message is used
			}
		}

		cilk_for (int j = 0; j < vertices.size(); j++) {
			for (int i = 0; i < nworkers; i++) {
				vertices[j] += localAccumulator[i][j];
			}
		}

		// Reset accumulator to be empty
		for (int w = 0; w < nworkers; w++) {
			std::fill(localAccumulator[w].begin(), localAccumulator[w].end(), 0);
		}
	}
}
