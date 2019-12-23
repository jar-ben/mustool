#ifndef BOOLEANSOLVER_H
#define BOOLEANSOLVER_H

#include "satSolvers/SatSolver.h"
#include <string>
#include <map>
#include <vector>
#include <unordered_map>

class BooleanSolver: public SatSolver{
public:
	BooleanSolver(std::string filename);
    	~BooleanSolver();
	std::vector<std::vector<int>> clauses;

        //model rotation
        int rotated_crits;
        std::vector<std::vector<int>> hitmap_pos;
        std::vector<std::vector<int>> hitmap_neg;
        std::vector<bool> flip_edges_computed;
        std::vector<std::vector<std::vector<int>>> flip_edges; // flip_edges[i][j][k] - i-th clause, j-th literal in the clause, k-th edge from i-th clause under j-th literal in the flip grap
        std::vector<std::vector<int>> flip_edges_flatten;
        void compute_flip_edges(int c);
	void criticals_rotation(std::vector<bool>& criticals, std::vector<bool> subset);
};
#endif
