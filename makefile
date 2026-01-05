### -- Host OS --
ifeq ($(OS),Windows_NT)
  ifeq ($(shell uname -s),) # not in a bash-like shell
	CLEANUP = del /F /Q
	MKDIR = mkdir
  else # in a bash-like shell, like msys
	CLEANUP = rm -f
	MKDIR = mkdir -p
  endif
	TARGET_EXTENSION=exe
else
	CLEANUP = rm -f
	MKDIR = mkdir -p
	TARGET_EXTENSION=out
endif

### -- Toolchain and flags --
SHELL = /bin/sh

# change extension depending on system (exe for windows, none for linux)
EXECUTABLE := mesosim

debug := 1
# XS_CFLAGS := -mcmodel=medium
# ?= means 'if not yet defined'
# variables used in definitions with '=' are expanded on use
BASE_CFLAGS = -std=c11
DEBUG_CFLAGS = -ggdb -O0 -Wall -Wextra -Wpedantic -Wshadow -Wconversion\
-DTEST -DUNITY_INCLUDE_CONFIG_H -I$(TEST_PATH)
OPT_CFLAGS = -O2
ALL_CFLAGS += $(CFLAGS) $(BASE_CFLAGS)

ifeq ($(debug), 1)
# add TEST_PATH to have access to unity_config.h 
	ALL_CFLAGS += $(DEBUG_CFLAGS)
else
	ALL_CFLAGS += $(OPT_CFLAGS)
endif
ALL_CFLAGS += -I. -I$(UNITY_PATH) -I$(SRC_PATH) -I$(INCLUDE_PATH)

# -v verbose
ALL_LDFLAGS += $(LDFLAGS) -lm

CC := gcc

# -c: only compile, don't link
COMPILE = $(CC) -c
LINK = $(CC)
DEPEND = $(CC) -MM -MG -MF
### -----

### -- Paths --
# := are expanded immediately
UNITY_PATH := unity/src
SRC_PATH := src
INCLUDE_PATH := include
TEST_PATH := test
BUILD_PATH := build
DEPEND_PATH := build/depends
OBJ_PATH := build/objs
RESULTS_PATH := build/results
# LIB_PATH := lib
BIN_PATH := bin

BUILD_PATHS := $(BUILD_PATH) $(DEPEND_PATH) $(OBJ_PATH) $(RESULTS_PATH)

TEST_SRC := $(wildcard $(TEST_PATH)/*.c)

$(info CFLAGS=$(ALL_CFLAGS))

# object files, use for mesosim target
# OBJS := $(addprefix $(BUILD_PATH)/,$(object_names))
# uses .c files in SRC_PATH to create list of object files in BUILD_PATH by pattern substitution
OBJS := $(patsubst $(SRC_PATH)/%.c,$(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))
SRC_MAIN := Mesosim.c
OBJS_MMAIN := $(filter-out $(OBJ_PATH)/$(subst .c,.o,$(SRC_MAIN)),$(OBJS))

# other headers to act as dependencies that don't have corresponding .c files
XS_HEADERS := State.h Defs.h Geometry.h
XS_HPATH := $(addprefix $(INCLUDE_PATH)/,$(XS_HEADERS))

### -- Test Results --
RESULTS := $(patsubst $(TEST_PATH)/Test%.c,$(RESULTS_PATH)/Test%.txt,$(TEST_SRC) )
PASSED = `grep :PASS $(RESULTS_PATH)/*.txt`
FAIL = `grep :FAIL $(RESULTS_PATH)/*.txt`
IGNORE = `grep :IGNORE $(RESULTS_PATH)/*.txt`
# SUMMARY = `tail -n2 $(RESULTS_PATH)/*.txt`
SUMMARY = `grep -T Tests $(RESULTS_PATH)/*.txt`

# $(info RESULTS is $(RESULTS))

.PHONY: test clean all
all: $(BIN_PATH)/$(EXECUTABLE).$(TARGET_EXTENSION)

# build binary using object files
# depends on object files and the headers that don't have corresponding objects
# requires paths as prerequisites so that folders get made if don't exist
$(BIN_PATH)/$(EXECUTABLE).$(TARGET_EXTENSION): $(BIN_PATH) $(BUILD_PATHS) $(OBJS) $(XS_HPATH)
	$(CC) $(ALL_CFLAGS) -o $@ $(OBJS) $(ALL_LDFLAGS)
	@touch print

# make object files based on .c files in SRC_PATH, looking for header files in INCLUDE_PATH
# wildcards needed on both sides, otherwise all .c files would be prerequisites
# $(BUILD_PATH)/%.o: $(SRC_PATH)/%.c $(INCLUDE_PATH)/%.h $(XS_HPATH)
# 	@mkdir -p $(BUILD_PATH) $(BIN_PATH)
# 	$(CC) $(ALL_CFLAGS) -I$(INCLUDE_PATH) -c -o $@ $<

# prints the filenames of files that were updated since the last 
# time the executable is built (tracked with the empty `print` file)
print: $(SRC_PATH)/* $(INCLUDE_PATH)/*
	@echo -ne " $(addsuffix \n,$?)"

test: $(BUILD_PATHS) $(RESULTS)
	@echo -e "-----------------------\nIGNORES:\n-----------------------"
	@echo "$(IGNORE)"
	@echo -e "-----------------------\nFAILURES:\n-----------------------"
	@echo "$(FAIL)"
	@echo -e "-----------------------\nPASSED:\n-----------------------"
	@echo "$(PASSED)"
	@echo -e "-----------------------\nSUMMARY:\n-----------------------"
	@echo "$(SUMMARY)"
	@echo -e "\nDONE"

$(RESULTS_PATH)/%.txt: $(BUILD_PATH)/%.$(TARGET_EXTENSION)
	./$< > $@ 2>&1

$(BUILD_PATH)/Test%.$(TARGET_EXTENSION): $(OBJS_MMAIN) $(OBJ_PATH)/Test%.o $(OBJ_PATH)/%.o $(OBJ_PATH)/unity.o #$(DEPEND_PATH)Test%.d
	$(LINK) -o $@ $^ $(ALL_LDFLAGS)

# :: (double-colon) rules are independent from each other
# TEST_PATH/TestSomething.c
$(OBJ_PATH)/%.o:: $(TEST_PATH)/%.c
	$(COMPILE) $(ALL_CFLAGS) $< -o $@

# SRC_PATH/Something.c
$(OBJ_PATH)/%.o:: $(SRC_PATH)/%.c $(INCLUDE_PATH)/%.h $(XS_HPATH)
	$(COMPILE) $(ALL_CFLAGS) $< -o $@

# UNITY_PATH/unity.c
$(OBJ_PATH)/%.o:: $(UNITY_PATH)/%.c $(UNITY_PATH)/%.h
	$(COMPILE) $(ALL_CFLAGS) $< -o $@

$(DEPEND_PATH)/%.d:: $(TEST_PATH)/%.c
	$(DEPEND) $@ $<

$(BUILD_PATHS) $(BIN_PATH):
	$(MKDIR) $@

# Clean build and bin directories
clean:
	$(CLEANUP) -r $(BUILD_PATH) $(BIN_PATH)
	$(CLEANUP) $(OBJ_PATH)/*.o
	$(CLEANUP) $(BUILD_PATH)/*.$(TARGET_EXTENSION)
	$(CLEANUP) $(RESULTS_PATH)/*.txt

.PRECIOUS: $(BUILD_PATH)/Test%.$(TARGET_EXTENSION)
.PRECIOUS: $(DEPEND_PATH)/%.d
.PRECIOUS: $(OBJ_PATH)/%.o
.PRECIOUS: $(RESULTS_PATH)/%.txt
