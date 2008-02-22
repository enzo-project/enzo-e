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
# ( ) dep
#
# Directories
#
# (*) data
# (*) io
#
#-----------------------------------------------------------------------

DIR = data io

all:   $(patsubst %,%-all,$(DIR))
clean: $(patsubst %,%-clean,$(DIR))
test:  $(patsubst %,%-test,$(DIR))
dep:   $(patsubst %,%-dep,$(DIR))
doc:
	doxygen

%-all: 
	$(MAKE) -C $*
%-clean: 
	-rm -f *~
	-rm -rf html latex
	$(MAKE) -C $* clean
%-test: 
	$(MAKE) -C $* test
%-dep: 
	$(MAKE) -C $* dep

