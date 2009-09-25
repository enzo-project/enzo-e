
TOP     = ../..

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: install-include
install-include:
	ln -sf `pwd`/test.hpp       $(INCLUDE)
	ln -sf `pwd`/test_unit.hpp  $(INCLUDE)

run-tests:
