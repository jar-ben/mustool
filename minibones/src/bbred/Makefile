GENCC = BBRed.cc
GENHH = BBRed.hh
ORIGCSRCS = $(filter-out $(GENCC), $(wildcard *.cc) ) \
            $(wildcard ../Lit*.cc) $(wildcard ../Read*.cc) ../minisat_aux.cc
CSRCS     = $(ORIGCSRCS) $(GENCC)
COBJS     = $(CSRCS:.cc=.o)

LIBD =
LIBS =
CXX?=g++

ifdef PROF
CFLAGS+=-pg
LNFLAGS+=-pg
CFLAGS+=-ggdb
endif

ifdef DBG
CFLAGS+=-O0
CFLAGS+=-ggdb
CFLAGS+=-DDBG
MSAT=libd
else
CFLAGS+=-DNDEBUG
CFLAGS+=-O3
MSAT=libr
endif

ifdef NOO
CFLAGS+=-O0
endif

CFLAGS += -Wall -DBUILDDATE='"$(BUILDDATE)"' -DDISTDATE='"$(DISTDATE)"'
CFLAGS += -DCHANGESET='"$(CHANGESET)"' -DRELDATE='"$(RELDATE)"'
CFLAGS+=-D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -Wno-parentheses -Wno-deprecated -D _MSC_VER
CFLAGS+=-I../minisat/ -I.. -I./generated
CFLAGS+=
CFLAGS+=-std=c++0x
LIBS+=-lminisat -lz 
LIBD+=-L../minisat/core

ifdef STATIC
CFLAGS+=-static
LNFLAGS+=-static
endif

NW_FILTER=-filter C++.filter

TARGET=bbred

.PHONY: m sources clean

all: sources $(TARGET)

$(TARGET):  m $(COBJS)
	@echo Linking: $@
	$(CXX)  -o $@ $(COBJS) $(LNFLAGS) $(LIBS) $(LIBD)

m:
	export MROOT=`pwd`/../minisat ; cd ../minisat/core; make CXX=$(CXX) LIB=minisat $(MSAT)

lib:  sources m $(COBJS)
	ar -rcsv lib$(TARGET).a `ls -1 *.o | grep -v 'main.o'`

## Build rule
%.o:	%.cc
	@echo Compiling: $@
	@$(CXX) $(CFLAGS) -c -o $@ $<


## Clean rule
clean:
	@rm -f  lib$(TARGET).a $(TARGET) $(TARGET).exe $(COBJS)
	@rm -f `cat .gitignore`
	@rm -f $(DEPFILE)
	@export MROOT=`pwd`/../minisat ; cd ../minisat/core; make CXX=$(CXX) clean

sources: $(GENCC) $(GENHH)

BBRed.cc: BBRed.hh
main.cc: BBRed.hh

BBRed.cc: BBRed.nw
	notangle -L -RBBRed.cc BBRed.nw >BBRed.cc

BBRed.hh: BBRed.nw
	notangle -L -RBBRed.hh BBRed.nw >BBRed.hh

doc: BBRed.pdf

%.pdf:  %.tex texheader.tex Makefile
	latexmk -pdf $<

%.tex:  %.nw
	noweave -x $(NW_FILTER) $< | \
	sed 's/documentclass{article}/documentclass{article}\\input{texheader.tex}/' >$@
#	noweave -x -index -autodefs c $(NW_FILTER) $< >$@  
