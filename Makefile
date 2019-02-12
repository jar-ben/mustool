DIR	= $(shell pwd)
CSRCS	= $(wildcard *.cpp)
COBJS	= $(CSRCS:.cpp=.o)

MINISAT	= $(DIR)/custom_minisat
MCSRCS	= $(wildcard $(MINISAT)/*.cc)
MCOBJS	= $(MCSRCS:.cc=.o)

LIBD 	= -L/usr/lib -L/usr/local/lib
LIBS 	= -lz -lspot -lz3
USR 	= /usr/include
INC 	= -I $(MINISAT) -I $(USR) -I /usr/local/include

CXX	= g++
#CFLAGS	= -O3 -w #-Wall
CFLAGS 	= -w -std=c++11

mvc: $(COBJS) $(MCOBJS)
	@echo Linking: $@
	$(CXX) -O3 -o $@ $(COBJS) $(MCOBJS) $(CFLAGS) $(INC) $(LIBD) $(LIBS) 

%.o: %.cpp
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

%.o: %.cc
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

print-%  : ; @echo $* = $($*)

clean:
	rm *.o
