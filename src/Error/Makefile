#
# ( ) liberror.a:
# ( ) clean:        Clean directory
# ( ) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

PACKAGE   = error
OBJ       = 

#=======================================================================
INC_LOCAL  = 
OBJ_EXTERN = 

OBJS = $(OBJ)
#=======================================================================

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: $(OBJS)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ lib*.a *.unit

#-----------------------------------------------------------------------

.PHONY: install

install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/error.hpp            $(INCLUDE)
	ln -sf `pwd`/error_exception.hpp  $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include




