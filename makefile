SHELL = /bin/sh

debug := 1
SRC_DIR := src
BUILD_DIR := build
INCLUDE_DIR := include
# LIB_DIR := lib
# TESTS_DIR := tests
BIN_DIR := bin

# change extension depending on system (exe for windows, none for linux)
EXECUTABLE := mesosim.exe

XS_CFLAGS := -mcmodel=medium
CFLAGS := -std=c11 -Wall -Wextra -Wpedantic $(XS_CFLAGS)
# -v verbose
LDFLAGS := -lm

ifeq ($(debug), 1)
	CFLAGS := $(CFLAGS) -ggdb -O0
else
	CFLAGS := $(CFLAGS) -O2
endif

# object files, use for mesosim target
# object_names := Mesosim.o Atoms.o FileIO.o Random.o Simulation_Aux.o Simulation.o Vector.o
# OBJS := $(addprefix $(BUILD_DIR)/,$(object_names))
OBJS := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o, $(wildcard $(SRC_DIR)/*.c))
XS_HEADERS := Defs.h Geometry.h
XS_HPATH := $(addprefix $(INCLUDE_DIR)/,$(XS_HEADERS))
# all: $(BIN_DIR)/$(EXECUTABLE)

# build binary using object files
# depends on object files and the headers that don't have corresponding objects
$(BIN_DIR)/$(EXECUTABLE): $(OBJS) $(XS_HPATH)
	gcc $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)
	@touch print

# make object files based on .c files in SRC_DIR, looking for header files in INCLUDE_DIR
# wildcards needed on both sides, otherwise all .c files would be prerequisites
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCLUDE_DIR)/%.h $(XS_HPATH)
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)
	gcc $(CFLAGS) -I$(INCLUDE_DIR) -o $@ -c $<

# Clean build and bin directories
clean:
	@rm -rf $(BUILD_DIR) $(BIN_DIR)

# prints the filenames of files that were updated since the last 
# time the executable is built (tracked with the empty `print` file)
print: $(SRC_DIR)/* $(INCLUDE_DIR)/*
	@echo -ne " $(addsuffix \n,$?)"
