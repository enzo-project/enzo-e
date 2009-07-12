#
# (x) clean:        Clean directory
# (^) dep:          Update dependencies
# (X) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

PACKAGE   = array
OBJ       = array_array_serial.o

#=======================================================================
OBJ_EXTERN = 

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: $(OBJS)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ out.make.dep Makefile.bak

#-----------------------------------------------------------------------

.PHONY: install
install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)   
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/array.hpp               $(INCLUDE)
	ln -sf `pwd`/array_array.hpp         $(INCLUDE)
	ln -sf `pwd`/array_array_serial.hpp  $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include

# DO NOT DELETE
