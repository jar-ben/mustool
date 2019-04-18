###################################################################################################
all:	r
-include config.mk

# Directory to store object files, libraries, executables, and dependencies:
BUILD_DIR      ?= build

# Target file names
MCSMUS      = mcsmus#       Name of Mcsmus main executable.

MCSMUS_CXXFLAGS = -I. -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -Wall -Wno-parentheses -Wextra -std=c++17
MCSMUS_LDFLAGS  = -Wall -lz 

SRCS = $(wildcard minisat/core/*.cc) $(wildcard minisat/simp/*.cc) $(wildcard minisat/utils/*.cc) \
	$(wildcard glucose/core/*.cc) $(wildcard glucose/simp/*.cc) $(wildcard glucose/utils/*.cc) \
	$(wildcard mcsmus/*.cc)
HDRS = $(wildcard minisat/mtl/*.h) $(wildcard minisat/core/*.h) $(wildcard minisat/simp/*.h) $(wildcard minisat/utils/*.h) \
	$(wildcard glucose/mtl/*.h) $(wildcard glucose/core/*.h) $(wildcard glucose/simp/*.h) $(wildcard glucose/utils/*.h) \
	$(wildcard mcsmus/*.hh)
OBJS = $(filter-out %Main.o, $(SRCS:.cc=.o))

r:	$(BUILD_DIR)/release/bin/$(MCSMUS)

## Build-type Compile-flags:
$(BUILD_DIR)/release/%.o:			MCSMUS_CXXFLAGS += -O3 -g

## Build-type Link-flags:
$(BUILD_DIR)/release/bin/$(MCSMUS):		MCSMUS_LDFLAGS += --static -g

## Executable dependencies
$(BUILD_DIR)/release/bin/$(MCSMUS):	 	$(foreach o,$(OBJS),$(BUILD_DIR)/release/$(o)) $(BUILD_DIR)/release/mcsmus/Main.o

## Compile rules
$(BUILD_DIR)/release/%.o:	%.cc
	@echo Compiling: $@
	mkdir -p $(dir $@)
	$(CXX) $(MCSMUS_CXXFLAGS) $(CXXFLAGS) -c -o $@ $< -MMD -MF $(BUILD_DIR)/release/$*.d

## Linking rule
$(BUILD_DIR)/release/bin/$(MCSMUS):
	@echo Linking Binary: $@
	mkdir -p $(dir $@)
	$(CXX) $^ $(MCSMUS_LDFLAGS) $(LDFLAGS) -o $@

clean:
	rm -f $(foreach t, release debug profile dynamic, $(foreach o, $(SRCS:.cc=.o), $(BUILD_DIR)/$t/$o)) \
          $(foreach t, release debug profile dynamic, $(foreach d, $(SRCS:.cc=.d), $(BUILD_DIR)/$t/$d)) \
	  $(foreach t, release debug profile dynamic, $(BUILD_DIR)/$t/bin/$(MCSMUS_CORE) $(BUILD_DIR)/$t/bin/$(MCSMUS)) \
