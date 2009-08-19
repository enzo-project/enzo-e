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
all: $(OBJS) tests

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ out.make.dep Makefile.bak *.unit
	-@rm parse parse.c *.out parse.tab.? lex.yy.c




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

parse: parse.y parse.l
	bison -d parse.y 
	flex parse.l
	gcc -o parse parse.tab.c lex.yy.c -lfl

tests: parse
	-parse < test.in         >& test.out
	-parse < test_list.in    >& test_list.out
	-parse < test_logical.in >& test_logical.out
	-parse < test_scalar.in  >& test_scalar.out
	-grep -i -l error *.out

#=======================================================================
#	INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include

# DO NOT DELETE
