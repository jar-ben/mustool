#include "SatSolver.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


class NuxmvHandle: public SatSolver{
public:
	int variables;
	std::vector<std::string> clauses_string;
	std::string model_file;
	std::string batch_file;
	std::string batch_prefix;
	int random_id;

	void build_nuxmv_model();
	void build_nuxmv_batch();
	void build_nuxmv_batch_file(std::vector<bool> &controls);

	std::string toString(std::vector<bool> &mus);
	NuxmvHandle(std::string filename);
	~NuxmvHandle();
	bool solve(std::vector<bool> &f, bool shrink = false, bool grow = false);
};
