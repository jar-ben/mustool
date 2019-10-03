DIR	= $(shell pwd)
MINISAT	= $(DIR)/custom_minisat
MCSMUS	= $(DIR)/mcsmus
MSAT	= libr

LIBD 	= -L/usr/lib -L/usr/local/lib
LIBS 	= -lz -lcryptominisat5
LIBS	+= -lstdc++fs
USR 	= /usr/include
INC 	= -I $(MCSMUS) -I $(MINISAT) -I $(USR) -I /usr/local/include -I $(DIR) -I $(MCSMUS) 

CSRCS	= $(wildcard *.cpp) $(wildcard $(DIR)/algorithms/*.cpp)
CSRCS	+= $(wildcard $(DIR)/satSolvers/*.cpp) $(wildcard $(DIR)/core/*.cpp) $(wildcard $(DIR)/counting/*.cpp)
COBJS	= $(CSRCS:.cpp=.o)

MCSRCS	= $(wildcard $(MINISAT)/*.cc)
MCOBJS	= $(MCSRCS:.cc=.o)

MCSMUS_SRCS = $(wildcard $(MCSMUS)/minisat/core/*.cc) $(wildcard $(MCSMUS)/minisat/simp/*.cc) $(wildcard $(MCSMUS)/minisat/utils/*.cc) \
		$(wildcard $(MCSMUS)/glucose/core/*.cc) $(wildcard $(MCSMUS)/glucose/simp/*.cc) $(wildcard $(MCSMUS)/glucose/utils/*.cc) \
		$(wildcard $(MCSMUS)/mcsmus/*.cc)
MCSMUS_OBJS = $(filter-out %Main.o, $(MCSMUS_SRCS:.cc=.o))

### 
# The following 3 variables control whether a support for individual constraint domains, SAT, SMT, LTL, should be build. 
USESAT = YES
USESMT = NO
USELTL = YES
###

USEMCSMUS = YES

CXX	= g++
CFLAGS 	= -w -std=c++17 -g
CFLAGS	+= -O3
CFLAGS	+= -D NDEBUG

ifeq ($(USEMCSMUS),YES)
	CFLAGS += -D UMCSMUS
else
	MCSMUS_OBJS = 
endif

ifeq ($(USESAT),NO)
	CFLAGS += -D NOSAT
endif
ifeq ($(USESMT),NO)
	CFLAGS += -D NOSMT
	CSRCS := $(filter-out $(DIR)/satSolvers/Z3Handle.cpp, $(CSRCS))
	COBJS := $(filter-out $(DIR)/satSolvers/Z3Handle.o, $(COBJS))
else
	LIBS   += -lz3
endif
ifeq ($(USELTL),NO)
	CFLAGS += -D NOLTL
	CSRCS := $(filter-out $(DIR)/satSolvers/SpotHandle.cpp, $(CSRCS))
	COBJS := $(filter-out $(DIR)/satSolvers/SpotHandle.o, $(COBJS))
	COBJS := $(filter-out $(DIR)/satSolvers/NuxmvHandle.o, $(COBJS))
	COBJS := $(filter-out $(DIR)/satSolvers/NuxmvHandle.cpp, $(COBJS))
else
	LIBS    += -lspot
endif

mvc: $(COBJS) $(MCOBJS) $(BONES_OBJS) $(MCSMUS_OBJS)
	@echo Linking: $@
	$(CXX) -o $@ $(COBJS) $(MCOBJS) $(MCSMUS_OBJS) $(CFLAGS) $(INC) $(LIBD) $(LIBS) 

%.o: %.cpp
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

%.o: %.cc
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) $(INC) -c -o $@ $<


print-%  : ; @echo $* = $($*)

clean:
	rm -f $(MCSMUS_OBJS)
	rm -f $(COBJS)
	rm -f $(MINISAT)/*.o

cleanCore:
	rm -f $(COBJS)
