
TOP     = ../../..
PACKAGE = user
OBJS    =  
TESTS   = 
EXES    = 
DIRS    = ppm ppml

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

 .PHONY: all
all: $(OBJS) $(TESTS) $(EXES)
	$(MAKE) -C ppm
	$(MAKE) -C ppml

#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	$(MAKE) -C ppm clean
	$(MAKE) -C ppml clean

