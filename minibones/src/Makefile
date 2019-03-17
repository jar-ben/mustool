DIR = $(shell pwd)
COBJS= $(DIR)/*.o
CXX?=g++

CFLAGS=-DNDEBUG
CFLAGS+=-O3


CFLAGS+=-I./minisat/
CFLAGS+=-I./picosat/
CFLAGS+=-std=c++0x
LIBS=-lz


LIBS+=-lminisat
LIBD=-L./minisat/core


all: $(COBJS)
	$(CXX) main.cc -o minibones $(COBJS) $(CFLAGS) $(LIBD) $(LIBS)

## Build rule
%.o:	%.cc
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) -c -o $@ $<

