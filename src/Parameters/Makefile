#
# (x) clean:           Clean directory
# (^) dep:             Update dependencies
# (X) install:         Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

PACKAGE   = parameters
OBJ       = parameters.o parse

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
	-@rm *.o *~ out.make.dep Makefile.bak *.unit
	-@rm parse parse.c


#-----------------------------------------------------------------------

.PHONY: install
install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)   
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/parameters.hpp               $(INCLUDE)

parse: parse.l
	flex -t parse.l > parse.c
	gcc -o parse parse.c -ll

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include

# DO NOT DELETE
