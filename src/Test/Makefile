# (X) install:         Install lib, include, and bin files
# (X) install-install: Install lib, include, and bin files
# 

TOP = ../..
SRC = ..
LIB = -L$(TOP)/lib -lcello

#=======================================================================
# TARGETS
#=======================================================================

all:

INCLUDE = $(TOP)/include

#=======================================================================

.PHONY: install
install: install-include $(OBJS)

.PHONY: install-include
install-include:
	ln -sf `pwd`/test.hpp       $(INCLUDE)
	ln -sf `pwd`/test_unit.hpp  $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include
