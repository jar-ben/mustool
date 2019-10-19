#include "core/misc.h"
#include "NuxmvHandle.h"
#include <algorithm>
#include <stdlib.h>

using namespace std;


NuxmvHandle::NuxmvHandle(string filename):SatSolver(filename){
	ifstream file;
	file.open(filename, ifstream::in);
	std::string tf;
	while (std::getline(file, tf)){
		clauses_string.push_back(tf);
	}
	dimension = clauses_string.size();

	//count number of variables
	stringstream cmd2;
	cmd2 << "grep -oP '(?<=p)[0-9]+' " << filename << " | awk 'NR == 1 { max=$1 } { if ($1>max) max=$1; } END {print max}'";
	string output = exec(cmd2.str().c_str());
	variables = std::stoi (output) + 1;
	random_id = random_number();

	build_nuxmv_model();
	build_nuxmv_batch();
}

NuxmvHandle::~NuxmvHandle(){
}

string NuxmvHandle::toString(vector<bool> &mus){
	stringstream result;
	for(int i = 0; i < dimension; i++){
		if(mus[i]){
			result << clauses_string[i] << "\n";

		}
	}
	return result.str();
}

void NuxmvHandle::build_nuxmv_model(){
	stringstream model, filen;
	model << "MODULE main" << std::endl << "VAR" << std::endl;
	for(int i = 0; i < variables; i++)
		model << "p" << i << ": boolean;" << std::endl;
	filen << "/tmp/model" << random_id << ".xmv";
	model_file = filen.str();
	std::ofstream file;
	file.open(model_file);
	file << model.str();
	file.close();	
}

void NuxmvHandle::build_nuxmv_batch(){
	stringstream batch, filen;
	batch << "read_model -i " << model_file << std::endl;
	batch << "go_msat" << std::endl; 
	batch << "set bmc_length 3" << std::endl;
	for(int i = 0; i < clauses_string.size(); i++)
		batch << "add_property -l -p \"" << clauses_string[i] << "\"" << std::endl;
	batch << "build_boolean_model" << std::endl;
	batch << "bmc_setup" << std::endl;
	batch_prefix = batch.str();
	filen << "/tmp/nuxmv_batch" << random_id;
	batch_file = filen.str();

	std::ofstream file;
	file.open(batch_file, ios_base::out);
	file << batch_prefix;
	file.close();		
}

void NuxmvHandle::build_nuxmv_batch_file(std::vector<bool> &controls){	
	std::stringstream batch_final;
	batch_final << "reqan_check_consistency -r \"";
	std::vector<string> ones;
	for(int i = 0; i < controls.size(); i++)
		if(controls[i])
			ones.push_back(std::to_string(i));
	if(!ones.empty())
		batch_final << ones[0];
	for(int i = 1; i < ones.size(); i++)
		batch_final << "," << ones[i];

	batch_final << "\" -e ic3" << std::endl << "quit";
		
	std::ofstream file;
	file.open(batch_file, ios_base::app);
	file << batch_final.str();
	file.close();		
}


bool NuxmvHandle::solve(std::vector<bool> &f, bool shrink, bool grow){
	if(std::count(f.begin(), f.end(), true) == 0)
		return true;
	checks++;
	bool ret = true;
	build_nuxmv_batch_file(f);
	std::stringstream output;
	output << "/tmp/nuxmv_result" << random_id;

	stringstream cmd;

	cmd << "../nuxmv/bin/nuXmv -int -source " << batch_file << " | tail -n 1 > " << output.str(); // << " > /dev/null";
	int status = system(cmd.str().c_str());


	std::ifstream ifs(output.str());
	std::string result( (std::istreambuf_iterator<char>(ifs) ),
		(std::istreambuf_iterator<char>()    ) );
	ret = (result.find("NOT consistent") == std::string::npos) && (result.find("consistent") != std::string::npos);

	//restore the batch file (remove last two lines)
	std::stringstream restore, temp;
	temp << "temp.txt" << random_id;
	restore << "head -n -2 " << batch_file << "> " <<temp.str() << "; cp " << temp.str() << " " << batch_file;
	system(restore.str().c_str());
	return ret;
}

