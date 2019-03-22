DIR	= $(shell pwd)
MINISAT	= $(DIR)/custom_minisat
BONES	= $(DIR)/minibones/src
MSAT	= libr

LIBD 	= -L/usr/lib -L/usr/local/lib
LIBD 	+= -L$(BONES)/minisat/build/release/lib
LIBS 	= -lz -lspot -lz3
LIBS	+= -lminisat
USR 	= /usr/include
INC 	= -I $(MINISAT) -I $(USR) -I /usr/local/include -I $(BONES)/minisat/ -I $(BONES)

CSRCS	= $(wildcard *.cpp)
COBJS	= $(CSRCS:.cpp=.o)

MCSRCS	= $(wildcard $(MINISAT)/*.cc)
MCOBJS	= $(MCSRCS:.cc=.o)

BCSRCS	= $(wildcard $(BONES)/*.cc)
BCOBJS	= $(BCSRCS:.cc=.o)

MINIBONESTYPE = MINIBONES_SRC

CXX	= g++
CFLAGS 	= -w -std=c++11 -g
CFLAGS	+= -O3
CFLAGS	+= -D NDEBUG
CFLAGS 	+= -D $(MINIBONESTYPE)


mvc: m $(COBJS) $(MCOBJS) $(BCOBJS)
	@echo Linking: $@
	$(CXX) -o $@ $(COBJS) $(MCOBJS) $(BCOBJS) $(CFLAGS) $(INC) $(LIBD) $(LIBS) 

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
	rm *.o
	rm $(MINISAT)/*.o
	rm $(BONES)/*.o
	find $(BONES)/minisat/ -name '*.o' | xargs rm -f
