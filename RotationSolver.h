#include <string>
#include <map>
#include <vector>
#include "custom_minisat/Solver.h"
#include <algorithm>


class RotationSolver{
public:
	int dimension; 
	int controls;
	std::vector<bool> controls_type;
	Minisat::Solver* solver;

	RotationSolver(int dim);
	~RotationSolver();

	void block_up(std::vector<bool> &f);
	void add_rotation_seed(std::vector<bool> &seed);
	void add_mus(
		std::vector<bool> &mus, 
		std::vector<std::vector<std::vector<int>>> &flip_edges, 
		std::vector<std::vector<int>> &parent_muses, 
		std::vector<std::vector<bool>> &muses,
		std::vector<std::vector<int>> muses_int,
		std::vector<bool> &critical);

	std::vector<bool> get_seed();
};
