DIR	= $(shell pwd)
MINISAT	= $(DIR)/custom_minisat
BONES	= $(DIR)/minibones/src
MCSMUS	= $(DIR)/mcsmus
MSAT	= libr

LIBD 	= -L/usr/lib -L/usr/local/lib
LIBD 	+= -L$(BONES)/minisat/build/release/lib
LIBS 	= -lz -lspot -lz3
LIBS	+= -lminisat -lstdc++fs
USR 	= /usr/include
INC 	= -I $(MCSMUS) -I $(MINISAT) -I $(USR) -I /usr/local/include -I $(BONES)/minisat/ -I $(BONES) -I $(DIR) -I $(MCSMUS) 

CSRCS	= $(wildcard *.cpp) $(wildcard $(DIR)/algorithms/*.cpp) $(wildcard $(DIR)/heuristics/*.cpp)
CSRCS	+= $(wildcard $(DIR)/satSolvers/*.cpp) $(wildcard $(DIR)/core/*.cpp)
COBJS	= $(CSRCS:.cpp=.o)

MCSRCS	= $(wildcard $(MINISAT)/*.cc)
MCOBJS	= $(MCSRCS:.cc=.o)

BCSRCS	= $(wildcard $(BONES)/*.cc)
BCOBJS	= $(BCSRCS:.cc=.o)

MCSMUS_SRCS = $(wildcard $(MCSMUS)/minisat/core/*.cc) $(wildcard $(MCSMUS)/minisat/simp/*.cc) $(wildcard $(MCSMUS)/minisat/utils/*.cc) \
		$(wildcard $(MCSMUS)/glucose/core/*.cc) $(wildcard $(MCSMUS)/glucose/simp/*.cc) $(wildcard $(MCSMUS)/glucose/utils/*.cc) \
		$(wildcard $(MCSMUS)/mcsmus/*.cc)
MCSMUS_OBJS = $(filter-out %Main.o, $(MCSMUS_SRCS:.cc=.o))

### 
# The following 3 variables control whether a support for individual constraint domains, SAT, SMT, LTL, should be build. 
USAT = YES
USMT = NO
ULTL = NO
###

CXX	= g++
CFLAGS 	= -w -std=c++17 -g
CFLAGS	+= -O3
CFLAGS	+= -D NDEBUG

ifeq ($(USAT),NO)
	CFLAGS += -D NOSAT
endif
ifeq ($(USMT),NO)
	CFLAGS += -D NOSMT
	CSRCS := $(filter-out $(DIR)/satSolvers/Z3Handle.cpp, $(CSRCS))
	COBJS := $(filter-out $(DIR)/satSolvers/Z3Handle.o, $(COBJS))
endif
ifeq ($(ULTL),NO)
	CFLAGS += -D NOLTL
	CSRCS := $(filter-out $(DIR)/satSolvers/SpotHandle.cpp, $(CSRCS))
	COBJS := $(filter-out $(DIR)/satSolvers/SpotHandle.o, $(COBJS))
endif

mvc: m $(COBJS) $(MCOBJS) $(BCOBJS) $(MCSMUS_OBJS)
	@echo Linking: $@
	$(CXX) -o $@ $(COBJS) $(MCOBJS) $(BCOBJS) $(MCSMUS_OBJS) $(CFLAGS) $(INC) $(LIBD) $(LIBS) 

%.o: %.cpp
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

%.o: %.cc
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

m:
	@echo Making Minisat
	export MROOT=$(BONES)/minisat ; cd $(BONES)/minisat; make CXX=$(CXX)

print-%  : ; @echo $* = $($*)

clean:
	rm -f $(MCSMUS_OBJS)
	rm -f $(COBJS)
	rm -f $(MINISAT)/*.o
	rm -f $(BONES)/*.o
	find $(BONES)/minisat/ -name '*.o' | xargs rm -f
