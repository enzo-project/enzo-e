#-----------------------------------------------------------------------
#
# @file    Makefile
# @brief   Top-level makefile
# @author  James Bordner
# @version 1.0
#
# Targets
#
# (*) all
# (*) clean
# (*) test
# (*) doc
# (*) dep
# (*) install
# (*) install-include
#
# Directories
#
# (*) Array
# (*) Disk
# (*) Error
# (*) Memory
# ( ) Method
# (*) Parameters
# (*) Performance
# (*) Test
#
#-----------------------------------------------------------------------

DIR = \
      Amr \
      Array \
      Control \
      Disk \
      Error \
      Memory \
      Method \
      Monitor \
      Parallel \
      Parameters \
      Performance \
      Test

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all %-all
all:       $(patsubst %,%-all,$(DIR))
%-all:  %-DEP
	-@($(MAKE)  -C $*)

#-----------------------------------------------------------------------

.PHONY: clean %-clean
clean:   $(patsubst %,%-clean,$(DIR))
%-clean:  %-DEP
	-@rm -f *~
	-@($(MAKE) -s -C $* clean)

#-----------------------------------------------------------------------

# Run tests in Test directory, and display results in Tests/*.unit

.PHONY: 
run-tests: all install $(patsubst %,%-run-tests,$(DIR))
	-@echo
	-@echo "=================================================="
	-@echo "FAILED TESTS:"
	-@echo "=================================================="
	-@echo
	-@grep FAIL */*unit
	-@echo

%-run-tests:
	$(MAKE) -s -C $* run-tests

#-----------------------------------------------------------------------

.PHONY: dep %-dep
dep:       $(patsubst %,%-dep,$(DIR))
%-dep: %-DEP
	$(MAKE) -s -C $* dep
.PHONY: %-DEP
%-DEP:
	-@(touch $*/DEPEND)

#-----------------------------------------------------------------------

.PHONY: install %-install
install:   $(patsubst %,%-install,$(DIR))
%-install: %-DEP
	-@($(MAKE) -s -C $* install)

#-----------------------------------------------------------------------

.PHONY: install-include %-install-include
install-include:   $(patsubst %,%-install-include,$(DIR))
%-install-include: %-DEP
	-@$(MAKE) -s -C $* install-include

#-----------------------------------------------------------------------

.PHONY: doc
doc:
	doxygen

#-----------------------------------------------------------------------

.PHONY: cccc
cccc:
	cccc */*.cpp

#-----------------------------------------------------------------------
# Package dependencies
#-----------------------------------------------------------------------

.PHONY: Amr
Amr: Disk
Amr: Error
Amr: Memory
Amr: Monitor

Array: Error
Array: Test

Control: Error
Control: Monitor
Control: Parameters
Control: Test

Disk: Array
Disk: Error
Disk: Test

Error: Test

Memory: Error
Memory: Test

Method: Disk
Method: Error
Method: Performance

Monitor: Performance

Parallel: Error
Parallel: Parameters
Parallel: Performance

Parameters: Error
Parameters: Test

Performance: Error
Performance: Memory
Performance: Test

