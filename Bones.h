#if !defined(ReMUS_Bones_h)
#define ReMUS_Bones_h

#include "minibones/src/ToolConfig.hh"
#include "minibones/src/CoreBased.hh"
#include "minibones/src/LitSet.hh"
#include <iostream>
#include "types_h.h"


class Bones{
public:
	ToolConfig config;
	Bones();
	~Bones();

	void set_options();	
	std::vector<int> get_backbones(std::vector<std::vector<int>> &cls);
	
};

#endif
