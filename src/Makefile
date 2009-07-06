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
# ( ) Performance
# (*) Test
#
#-----------------------------------------------------------------------

DIR = Array Disk Test Performance Error

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
	-@rm -rf html latex
	$(MAKE) -C $* clean

#-----------------------------------------------------------------------

.PHONY: test %-test
test:      $(patsubst %,%-test,$(DIR))
%-test: %-DEP
	$(MAKE) -C $* test

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
	cp doxygen.css html

#-----------------------------------------------------------------------

.PHONY: cccc
cccc:
	cccc */*.cpp
