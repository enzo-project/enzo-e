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
# (*) run-tests
# (*) dep
# (*) install
# (*) install-include
# (*) doc
# (*) cccc
#
#-----------------------------------------------------------------------

DIR = \
      Amr \
      Array \
      Disk \
      Distribute \
      Error \
      Field \
      Memory \
      Method \
      Monitor \
      Parallel \
      Parameters \
      Particles \
      Performance \
      Schedule \
      Simulation \
      Task \
      Test

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all %-component
all: install-include dep $(patsubst %,%-component,$(DIR))
%-component:  %-DEP
	-@($(MAKE)  -C $*)

#-----------------------------------------------------------------------

.PHONY: clean %-clean
clean:   $(patsubst %,%-clean,$(DIR))
%-clean:  %-DEP
	-@rm -f *~
	-@($(MAKE) -s -C $* clean)

#-----------------------------------------------------------------------

# Run tests in Test directory, and display results in Tests/*.unit

.PHONY: run-tests
run-tests: all install $(patsubst %,%-run-tests,$(DIR))
	-@echo
	-@echo "=================================================="
	-@echo "FAILED TESTS:"
	-@echo "=================================================="
	-@echo
	-@$(MAKE) -s failed-tests

%-run-tests:
	$(MAKE) -s -C $* run-tests

#-----------------------------------------------------------------------

.PHONY: failed-tests
failed-tests:
	-@grep FAIL `find . -name "*.unit"` | awk '{print $$2,$$3}'
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


Amr-component: Disk-component
Amr-component: Error-component
Amr-component: Memory-component
Amr-component: Monitor-component

Array-component: Error-component
Array-component: Test-component

Schedule-component: Error-component
Schedule-component: Monitor-component
Schedule-component: Parameters-component
Schedule-component: Simulation-component
Schedule-component: Test-component

Disk-component: Array-component
Disk-component: Error-component
Disk-component: Test-component

Error-component: Test-component

Memory-component: Error-component
Memory-component: Performance-component
Memory-component: Test-component

Method-component: Disk-component
Method-component: Error-component
Method-component: Performance-component
Method-component: Monitor-component

Monitor-component: Performance-component
Monitor-component: Test-component

Parallel-component: Error-component
Parallel-component: Parameters-component
Parallel-component: Performance-component

Parameters-component: Error-component
Parameters-component: Test-component

Performance-component: Error-component
Performance-component: Memory-component
Performance-component: Test-component

Simulation-component: Amr-component
Simulation-component: Field-component
Simulation-component: Method-component
Simulation-component: Parallel-component
Simulation-component: Parameters-component


