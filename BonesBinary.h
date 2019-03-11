#if defined(MINIBONES_BIN) && !defined(ReMUS_BonesBinary_h)
#define ReMUS_BonesBinary_h

#include <iostream> 
#include <string> 
#include <sstream> 

class Bones{
public:
	Bones();
	~Bones();
	std::vector<int> get_backbones(std::vector<std::vector<int>> &cls);
};

Bones::Bones(){}
Bones::~Bones(){}

std::vector<int> Bones::get_backbones(std::vector<std::vector<int>> &cls){
	stringstream ss;
	int clauses = cls.size();
	int vars = 0;;
        for(auto &cl: cls){
                for(auto l: cl){
                        ss << l << " ";
			if(l >= vars || -1 * l >= vars) vars = l;
		}
                ss << "0" << "\n";
        }

        ofstream file;
        stringstream outputf;
  	int salt = 1;//TODO: change it to something random
  	outputf << "minibones_export" << salt << ".cnf";
        string filename = outputf.str();
        file.open(filename);
	file << "p cnf " << vars << " " << clauses << "\n";
        file << ss.str();
        file.close();

        //compute the backbone (invoke minibones)
        stringstream cmd;
        cmd << "./minibones_binary -e -i -c 100 " << filename << " 2> /dev/null";
        string output = exec(cmd.str().c_str());
        stringstream oss(output);
        string line;
	vector<int> backbone_lits;
	while (getline(oss, line, '\n'))
        {
                if (line[0] != 'v')
                        continue;
                istringstream is(line);
                is.ignore(1,' '); //ignore the initial 'v'
                int lit;
                while(is >> lit){
                        if(lit == 0) continue;
			backbone_lits.push_back(lit);
                }
        }
        remove(filename.c_str());
        return backbone_lits;
}

	
#endif
