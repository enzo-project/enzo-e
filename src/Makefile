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

DIR = Array Disk Error Memory Parameters Performance Error Test
# Method

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all %-all
all:       $(patsubst %,%-all,$(DIR))
%-all:  %-DEP
	$(MAKE) -C $*

#-----------------------------------------------------------------------

.PHONY: clean %-clean
clean:     $(patsubst %,%-clean,$(DIR))
%-clean:  %-DEP
	-@rm -f *~
	$(MAKE) -C $* clean

#-----------------------------------------------------------------------

# Run tests in Test directory, and display results in Tests/*.unit

.PHONY: 
test:   all install
	$(MAKE) -C Test test
	$(MAKE) -C Memory test
	-@cat Test/*.unit
	-@cat Memory/*.unit

#-----------------------------------------------------------------------

.PHONY: dep %-dep
dep:       $(patsubst %,%-dep,$(DIR))
%-dep: %-DEP
	$(MAKE) -C $* dep
.PHONY: %-DEP
%-DEP:
	touch $*/DEPEND

#-----------------------------------------------------------------------

.PHONY: install %-install
install:   $(patsubst %,%-install,$(DIR))
%-install: %-DEP
	$(MAKE) -C $* install

#-----------------------------------------------------------------------

.PHONY: install-include %-install-include
install-include:   $(patsubst %,%-install-include,$(DIR))
%-install-include: %-DEP
	$(MAKE) -C $* install-include

#-----------------------------------------------------------------------

.PHONY: doc
doc:
	doxygen

#-----------------------------------------------------------------------

.PHONY: cccc
cccc:
	cccc */*.cpp
