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
#
# Directories
#
# (*) array
# (*) test
# (*) float
# (*) io
#
#-----------------------------------------------------------------------

DIR = array test float io

all:   $(patsubst %,%-all,$(DIR))
clean: $(patsubst %,%-clean,$(DIR))
test:  $(patsubst %,%-test,$(DIR))

%-all: 
	$(MAKE) -C $*
%-clean: 
	-rm -f *~
	-rm -rf html latex
	$(MAKE) -C $* clean
%-test: 
	$(MAKE) -C $* test
doc:
	doxygen

