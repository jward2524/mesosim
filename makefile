# Makefile for Mesosim project
# release - optimized build
# debug - debug build with symbols and no optimization
# test - build with unity tests
# dbtest - debug build with unity tests

SHELL := /bin/sh
# http://redsymbol.net/articles/unofficial-bash-strict-mode/
.SHELLFLAGS := -eu -c
.DELETE_ON_ERROR:
# MAKEFLAGS += --warn-undefined-variables # this file uses a lot of undefined variables
MAKEFLAGS += --no-builtin-rules

# := are expanded immediately
# ?= means 'if not yet defined'
# variables used in definitions with '=' are expanded on use

# --- Host OS ---
ifeq ($(OS),Windows_NT)
    ifeq ($(shell uname -s),) # not in a bash-like shell
		cleanup = del /F /Q
		mkdir = mkdir
    else # in a bash-like shell, like msys
		cleanup = rm -f
		mkdir = mkdir -p
    endif
	TARGET_EXTENSION = exe
else
	cleanup = rm -f
	mkdir = mkdir -p
	TARGET_EXTENSION = out
endif

# gcc - : read from stdin
# no promises that core dumping will work if this test passes
# gcore may still fail due to security settings (e.g. yama ptrace_scope)
DUMP_CORE := $(and \
	$(shell echo 'int main(){fork();}' | \
	$(CC) -x c -o a.$(EXECUTABLE) - >/dev/null 2>&1 \
	&& rm a.$(EXECUTABLE) && echo 1), \
	$(shell gdb -batch-silent && echo 1)\
)
DUMPFLAG := $(if $(DUMP_CORE),-DDUMP_CORE)
# echo 'int main(){fork();}' | gcc -x c -o a.exe - >/dev/null 2>&1 && rm a.exe && echo yes

# --- Toolchain ---
CC := gcc

# -c: only compile, don't link
COMPILE = $(CC) -c
LINK = $(CC)
DEPEND = $(CC) -MM -MG -MF

TYPES := release debug test dbtest docs
DEFAULT_TYPE := release

# Deduce TYPE from a file path like build/<TYPE>/..., else empty string
ifndef TYPE
    ifeq ($(filter $(MAKECMDGOALS),$(TYPES)),)
		TYPE ?= $(firstword $(foreach t,$(TYPES), \
				$(if $(findstring build/$(t)/,$(MAKECMDGOALS)),$(t)) ))
#         $(info Deduced TYPE=$(TYPE) from MAKECMDGOALS=$(MAKECMDGOALS))
    endif
endif

TYPE ?= $(DEFAULT_TYPE)
DEBUG ?= $(if $(filter $(TYPE),debug dbtest),1,0)

# --- Flags ---
BASE_CFLAGS ?= -std=c11
# add test_path to find unity_config.h
# needed by unity_internals.h
unity_CFLAGS = -DTEST -DUNITY_INCLUDE_CONFIG_H -I$(test_path) -I$(unity_path)

ifeq ($(DEBUG),1)
	WARNFLAGS ?= -Wall -Wextra -Wpedantic -Wshadow -Wconversion
	OPTFLAGS ?= -O0
	DEBUGFLAGS ?= -ggdb
else
	BASE_CFLAGS += -DNDEBUG
	OPTFLAGS ?= -O2
	WARNFLAGS ?=
endif

debug_CFLAGS := $(DEBUGFLAGS) $(OPTFLAGS) $(WARNFLAGS)
release_CFLAGS := $(OPTFLAGS) $(WARNFLAGS)
test_CFLAGS = ${release_CFLAGS} ${unity_CFLAGS}
dbtest_CFLAGS = ${debug_CFLAGS} ${unity_CFLAGS}

INCLUDE_CFLAGS = -I$(src_path) -I$(include_path)

# CFLAGS can be overridden from the command line
ALL_CFLAGS = \
	$(BASE_CFLAGS) \
	$($(TYPE)_CFLAGS) \
	$(DUMPFLAG) \
	$(INCLUDE_CFLAGS) \
	$(CFLAGS)

# -v verbose
ALL_LDFLAGS = -lm $(LDFLAGS)

# --- Paths ---
# extension changes depending on system (exe for windows, .out for linux)
EXECUTABLE := mesosim

unity_path := unity/src
src_path := src
include_path := include
test_path := test
build_path := build
docs_path := docs
# depend_path := build/depends
# # LIB_PATH := lib

build_paths := $(build_path) $(depend_path)

# --- Build type paths ---
obj_path := $(build_path)/$(TYPE)/objs
bin_path := $(build_path)/$(TYPE)

ifneq (,$(filter $(TYPE),test dbtest))
	results_path := $(build_path)/$(TYPE)/results
endif
build_paths += $(bin_path) $(obj_path) $(results_path)

# --- Source and object files ---
# uses .c files in src_path to create list of object files in build_path by pattern substitution
# evaluated on use, after TYPE is set
objs = $(patsubst $(src_path)/%.c,$(obj_path)/%.o,$(wildcard $(src_path)/*.c))

src_main := Mesosim.c
objs_mmain = $(filter-out $(obj_path)/$(subst .c,.o,$(src_main)),$(objs))

# other headers to act as dependencies that don't have corresponding .c files
xs_headers := State.h Defs.h Geometry.h
xs_hpath := $(addprefix $(include_path)/,$(xs_headers))

test_src := $(wildcard $(test_path)/*.c)

# --- Build targets ---
# $(info CFLAGS=$(ALL_CFLAGS))
.PHONY: $(TYPES) clean all
.DEFAULT_GOAL := release
all: $(TYPES)

# unset variables are empty strings
do_build_deps = $(bin_path)/$(EXECUTABLE).$(TARGET_EXTENSION)
ifneq (,$(filter $(TYPE),test))
	do_build_deps = do-test
endif
ifneq (,$(filter $(TYPE),dbtest))
	do_build_deps = $(patsubst $(test_path)/Test%.c,$(bin_path)/Test%.$(TARGET_EXTENSION),$(test_src) )
endif
ifneq (,$(filter $(TYPE),docs))
	do_build_deps = do-docs
endif

$(TYPES):
	@echo -e "Building type: $@\n"
	@"$(MAKE)" TYPE=$@ do-build

do-build: $(do_build_deps) | $(build_paths)
# 	$(info TYPE=$(TYPE))
	@touch print-$(TYPE)

do-docs: $(bin_path)/InputMrk.$(TARGET_EXTENSION) | $(build_paths)
	@./$< > docs/Commands.md
	@echo "Documentation generated at docs/Commands.md"

# --- Test results ---
results = $(patsubst $(test_path)/Test%.c,$(results_path)/Test%.txt,$(test_src) )
passed = `grep :PASS $(results_path)/*.txt`
fail = `grep :FAIL $(results_path)/*.txt`
ignore = `grep :IGNORE $(results_path)/*.txt`
# summary = `tail -n2 $(results_path)/*.txt`
summary = `grep -T Tests $(results_path)/*.txt | column -t`

do-test: $(results) | $(build_paths)
# 	@echo -e "-----------------------\nPASSED:\n-----------------------"
# 	@echo "$(passed)"
	@echo -e "-----------------------\nIGNORES:\n-----------------------"
	@echo "$(ignore)"
	@echo -e "-----------------------\nFAILURES:\n-----------------------"
	@echo "$(fail)"
	@echo -e "-----------------------\nSUMMARY:\n-----------------------"
	@echo "$(summary)"
	@echo -e "\nDONE"

# --- Build rules ---
# build binary using object files
# depends on object files and the headers that don't have corresponding objects
$(bin_path)/$(EXECUTABLE).$(TARGET_EXTENSION): $(objs) $(xs_hpath) | $(bin_path)
	$(CC) $(ALL_CFLAGS) -o $@ $(objs) $(ALL_LDFLAGS)

$(results_path)/%.txt: $(bin_path)/%.$(TARGET_EXTENSION) | $(results_path)
	./$< > $@ 2>&1

$(bin_path)/Test%.$(TARGET_EXTENSION): $(objs_mmain) $(obj_path)/Test%.o $(obj_path)/%.o $(obj_path)/unity.o $(obj_path)/TUtils.o | $(bin_path)
	$(LINK) -o $@ $^ $(ALL_LDFLAGS)

# :: (double-colon) rules are independent from each other
# test_path/TestSomething.c
$(obj_path)/%.o:: $(test_path)/%.c $(test_path)/TUtils.h $(xs_hpath) | $(obj_path)
	$(COMPILE) $(ALL_CFLAGS) $< -o $@

# src_path/Something.c
$(obj_path)/%.o:: $(src_path)/%.c $(include_path)/%.h $(xs_hpath) | $(obj_path)
	$(COMPILE) $(ALL_CFLAGS) $< -o $@

# unity_path/unity.c
$(obj_path)/%.o:: $(unity_path)/%.c $(unity_path)/%.h $(test_path)/unity_config.h | $(obj_path)
	$(COMPILE) $(ALL_CFLAGS) $< -o $@

# $(depend_path)/%.d: $(test_path)/%.c
# 	$(DEPEND) $(ALL_CFLAGS) $@ $<

$(build_paths):
	$(mkdir) $@

$(bin_path)/InputMrk.$(TARGET_EXTENSION): $(obj_path)/InputMrk.o $(objs_mmain) | $(bin_path)
	$(CC) $(ALL_CFLAGS) $(objs_mmain) $< -o $@ $(ALL_LDFLAGS)

$(obj_path)/InputMrk.o: $(docs_path)/InputMrk.c | $(obj_path)
	$(COMPILE) $(ALL_CFLAGS) $< -o $@

# --- Print updated files ---
# prints the filenames of files that were updated since the last 
# time the executable is built (tracked with the empty `print-TYPE` file)
print-%: $(src_path)/* $(include_path)/* $(test_path)/*
	@echo -ne " $(addsuffix \n,$?)"

.PHONY: print-all print
# prints the updated files for the most recent build type
print:
	@file=$$(ls -tp print* | grep -v / | head -n 1); \
	type=$${file##print-}; \
	echo "Most recent type: $$type"; \
	"$(MAKE)" -s "$$file"

# prints the updated files for all build types
print-all:
	@for t in $(TYPES); do \
		echo "=== $$t ==="; \
		"$(MAKE)" -s print-$$t; \
	done

# --- Cleaning ---
# Clean build and bin directories
clean:
	$(cleanup) -r $(build_path)
	$(cleanup) print-*
# 	$(cleanup) $(obj_path)/*.o
# 	$(cleanup) $(build_path)/*.$(TARGET_EXTENSION)
# 	$(cleanup) $(results_path)/*.txt

# prevent deletion of intermediate files
.PRECIOUS: $(bin_path)/Test%.$(TARGET_EXTENSION)
.PRECIOUS: $(depend_path)/%.d
.PRECIOUS: $(obj_path)/%.o
.PRECIOUS: $(results_path)/%.txt
