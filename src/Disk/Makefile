#
# (x) clean:        Clean directory
# (^) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

PACKAGE   = disk
OBJ       = hdf5.o

#=======================================================================
OBJ_EXTERN = $(patsubst %,$(SRC)/%,$(ARRAY_OBJS))

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: $(OBJS)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ out.make.dep Makefile.bak lib*.a *.hdf5

#-----------------------------------------------------------------------

.PHONY: install
install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/disk.hpp       $(INCLUDE)
	ln -sf `pwd`/disk_hdf5.hpp $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include
include $(SRC)/MAKE/Make.objects



