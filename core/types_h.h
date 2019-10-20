#ifndef TYPES_H
#define TYPES_H

#include <vector>

typedef std::vector<bool> Formula;
typedef std::vector<int> Block;

class MUS{
	public:
	MUS(Formula& f, float d, int idc, int seed_d = 1){
		bool_mus = f;
		for(int i = 0; i < f.size(); i++)
			if(f[i])
				int_mus.push_back(i);
		without_crits = int_mus;		
		dimension = int_mus.size();
		duration = d;
		id = idc;
		seed_dimension = seed_d;
	}
	Formula bool_mus;
	int seed_dimension;
	std::vector<int> int_mus;
	std::vector<int> without_crits;
	int dimension;
	float duration;
	int id;
};


#endif
