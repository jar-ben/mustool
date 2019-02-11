DIR = $(shell pwd)
LIB = -L /usr/lib -L /usr/local/lib -L $(DIR)/include
LIB += -L $(DIR)/minibones/src/minisat/core
LINK = -lz -lspot
MINISAT = $(DIR)/custom_minisat
SOURCE = $(DIR)/*.cpp
#SOURCE += $(MINISAT)/*.cc
USR = /usr/include
INC = -I $(USR) -I /usr/local/include
#INC += -I $(MINISAT)
COBJS= $(DIR)/minibones/src/*.o
CXX?=g++

all:
	$(CXX) -O3 -g $(INC) $(LIB) $(LINK) $(SOURCE) -w -std=c++11 -o mvc $(COBJS) -lz3 -lspot -lminisat
