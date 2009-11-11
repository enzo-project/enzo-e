
TOP     = ../..
PACKAGE = amr
OBJS    = 
TESTS   = 
DIRS    = tree

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

# override all target to compile $(DIRS)

 .PHONY: all
all: $(OBJS) $(TESTS) $(EXES)
	$(MAKE) -C tree

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/amr.hpp       $(INCLUDE)
	$(MAKE) -C tree install-include

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	-@$(MAKE) -C tree clean

