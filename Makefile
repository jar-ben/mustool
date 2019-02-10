DIR = $(shell pwd)
LIB = -L /usr/lib -L /usr/local/lib -L $(DIR)/include
LINK = -lz -lspot -lmuser2_api
MINISAT = $(DIR)/custom_minisat
SOURCE = $(DIR)/*.cpp $(MINISAT)/*.cc
USR = /usr/include
INC = -I $(MINISAT) -I $(USR) -I /usr/local/include

all:
	g++ -O3 -g $(INC) $(LIB) $(LINK) $(SOURCE) -w -std=c++11 -o mvc -lz3 -lspot -lmuser2_api
