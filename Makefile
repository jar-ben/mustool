DIR	= $(shell pwd)
MINISAT	= $(DIR)/custom_minisat
BONES	= $(DIR)/minibones/src

LIBD 	= -L/usr/lib -L/usr/local/lib
LIBD 	+= -L$(BONES)/minisat/core
LIBS 	= -lz -lspot -lz3
LIBS	+= -lminisat
USR 	= /usr/include
INC 	= -I $(MINISAT) -I $(USR) -I /usr/local/include -I $(BONES)/minisat/ -I $(BONES)

CSRCS	= $(wildcard *.cpp)
COBJS	= $(CSRCS:.cpp=.o)

MCSRCS	= $(wildcard $(MINISAT)/*.cc)
MCOBJS	= $(MCSRCS:.cc=.o)

BCOBJS	= $(wildcard $(BONES)/*.o)


CXX	= g++
CFLAGS 	= -w -std=c++11 -g
CFLAGS	+= -O3
CFLAGS	+= -D NDEBUG


mvc: $(COBJS) $(MCOBJS) $(BCOBJS)
	@echo Linking: $@
	$(CXX) -o $@ $(COBJS) $(MCOBJS) $(BCOBJS) $(CFLAGS) $(INC) $(LIBD) $(LIBS) 

%.o: %.cpp
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

%.o: %.cc
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

print-%  : ; @echo $* = $($*)

clean:
	rm *.o
