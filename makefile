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
CFLAGS ?= -std=c11
ifeq ($(debug), 1)
	CFLAGS += -ggdb -O0 -Wall -Wextra -Wpedantic -DTEST
else
	CFLAGS ?= -O2
endif
CFLAGS += -I. -I$(UNITY_PATH) -I$(SRC_PATH) -I$(INCLUDE_PATH)

# -v verbose
LDFLAGS += -lm

CC := gcc

# -c: only compile, don't link
COMPILE = $(CC) -c
LINK = $(CC)
DEPEND = $(CC) -MM -MG -MF
### -----

### -- Paths --
UNITY_PATH = unity/src
SRC_PATH = src
INCLUDE_PATH = include
TEST_PATH = test
BUILD_PATH = build
DEPEND_PATH = build/depends
OBJ_PATH = build/objs
RESULTS_PATH = build/results
# LIB_PATH = lib
BIN_PATH = bin

BUILD_PATHS = $(BUILD_PATH) $(DEPEND_PATH) $(OBJ_PATH) $(RESULTS_PATH)

TEST_SRC = $(wildcard $(TEST_PATH)/*.c)

RESULTS = $(patsubst $(TEST_PATH)/Test%.c,$(RESULTS_PATH)/Test%.txt,$(TEST_SRC) )
PASSED = `grep -s PASS $(RESULTS_PATH)*.txt`
FAIL = `grep -s FAIL $(RESULTS_PATH)*.txt`
IGNORE = `grep -s IGNORE $(RESULTS_PATH)*.txt`

# object files, use for mesosim target
# OBJS := $(addprefix $(BUILD_PATH)/,$(object_names))
# uses .c files in SRC_PATH to create list of object files in BUILD_PATH by pattern substitution
OBJS := $(patsubst $(SRC_PATH)/%.c,$(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))

# other headers to act as dependencies that don't have corresponding .c files
XS_HEADERS := Common.h Defs.h Geometry.h
XS_HPATH := $(addprefix $(INCLUDE_PATH)/,$(XS_HEADERS))
### -----

all: $(BIN_PATH)/$(EXECUTABLE).$(TARGET_EXTENSION)

# build binary using object files
# depends on object files and the headers that don't have corresponding objects
# requires paths as prerequisites so that folders get made if don't exist
$(BIN_PATH)/$(EXECUTABLE).$(TARGET_EXTENSION): $(BIN_PATH) $(BUILD_PATHS) $(OBJS) $(XS_HPATH)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)
	@touch print

# make object files based on .c files in SRC_PATH, looking for header files in INCLUDE_PATH
# wildcards needed on both sides, otherwise all .c files would be prerequisites
# $(BUILD_PATH)/%.o: $(SRC_PATH)/%.c $(INCLUDE_PATH)/%.h $(XS_HPATH)
# 	@mkdir -p $(BUILD_PATH) $(BIN_PATH)
# 	$(CC) $(CFLAGS) -I$(INCLUDE_PATH) -c -o $@ $<

# prints the filenames of files that were updated since the last 
# time the executable is built (tracked with the empty `print` file)
print: $(SRC_PATH)/* $(INCLUDE_PATH)/*
	@echo -ne " $(addsuffix \n,$?)"

.PHONY: test
test: $(BUILD_PATHS) $(RESULTS)
	@echo "-----------------------\nIGNORES:\n-----------------------"
	@echo "$(IGNORE)"
	@echo "-----------------------\nFAILURES:\n-----------------------"
	@echo "$(FAIL)"
	@echo "-----------------------\nPASSED:\n-----------------------"
	@echo "$(PASSED)"
	@echo "\nDONE"

$(RESULTS_PATH)/%.txt: $(BUILD_PATH)/%.$(TARGET_EXTENSION)
	-./$< > $@ 2>&1

$(BUILD_PATH)/Test%.$(TARGET_EXTENSION): $(OBJ_PATH)/Test%.o $(OBJ_PATH)/%.o $(OBJ_PATH)/unity.o #$(DEPEND_PATH)Test%.d
	$(LINK) -o $@ $^

# :: (double-colon) rules are independent
# TEST_PATH/TestSomething.c
$(OBJ_PATH)/%.o:: $(TEST_PATH)/%.c $(INCLUDE_PATH)/%.h $(XS_HPATH)
	$(COMPILE) $(CFLAGS) $< -o $@

# SRC_PATH/Something.c
$(OBJ_PATH)/%.o:: $(SRC_PATH)/%.c $(INCLUDE_PATH)/%.h $(XS_HPATH)
	$(COMPILE) $(CFLAGS) $< -o $@

# UNITY_PATH/unity.c
$(OBJ_PATH)/%.o:: $(UNITY_PATH)/%.c $(UNITY_PATH)/%.h
	$(COMPILE) $(CFLAGS) $< -o $@

$(DEPEND_PATH)/%.d:: $(TEST_PATH)/%.c
	$(DEPEND) $@ $<

$(BUILD_PATHS) $(BIN_PATH):
	$(MKDIR) $@

# $(DEPEND_PATH):
# 	$(MKDIR) $(DEPEND_PATH)

# $(OBJ_PATH):
# 	$(MKDIR) $(OBJ_PATH)

# $(RESULTS_PATH):
# 	$(MKDIR) $(RESULTS_PATH)

# Clean build and bin directories
.PHONY: clean
clean:
	$(CLEANUP) -r $(BUILD_PATH) $(BIN_PATH)
	$(CLEANUP) $(OBJ_PATH)/*.o
	$(CLEANUP) $(BUILD_PATH)/*.$(TARGET_EXTENSION)
	$(CLEANUP) $(RESULTS_PATH)/*.txt

.PRECIOUS: $(BUILD_PATH)/Test%.$(TARGET_EXTENSION)
.PRECIOUS: $(DEPEND_PATH)/%.d
.PRECIOUS: $(OBJ_PATH)/%.o
.PRECIOUS: $(RESULTS_PATH)/%.txt
