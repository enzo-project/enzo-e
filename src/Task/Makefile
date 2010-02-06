
TOP     = ../..
PACKAGE = task
OBJS    = task.o
TESTS   = 
DIRS    = 

include $(TOP)/src/MAKE/Make.include

#=======================================================================

# override all target to compile $(DIRS)

 .PHONY: all
all: $(OBJS) $(TESTS) $(EXES)

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/task.hpp               $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:

