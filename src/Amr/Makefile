
TOP     = ../..
PACKAGE = amr
OBJS    = 
TESTS   = 
DIRS    = tree

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

# override all target to compile ppm and ppml
 .PHONY: all
all: $(OBJS) $(TESTS) $(EXES)
	$(MAKE) -C tree

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	-@$(MAKE) -C tree clean

