#
# (x) clean:           Clean directory
# (^) dep:             Update dependencies
# (X) install:         Install lib, include, and bin files
# 

TOP = ../..
SRC = ..
LIB = -L$(TOP)/lib -lcello

PACKAGE = parameters
OBJS    = parameters.o
TESTS   = test_parameters test_parse

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: $(OBJS) $(TESTS)

#-----------------------------------------------------------------------

.PHONY: run-tests
run-tests: $(TESTS)
	-./test_parameters
	-./test_parse < test_list.in    > test_list.out
	-./test_parse < test_logical.in > test_logical.out
	-./test_parse < test_scalar.in  > test_scalar.out
	-./test_parse < shock-pool-2.0.in > shock-pool-2.0.out
	-./test_parse < implosion-1.0.in > implosion-1.0.out
	-grep -i -l error *.out

#-----------------------------------------------------------------------

test_parameters: test_parameters.o
	$(CXX) -o test_parameters test_parameters.o $(LIB)

test_parse: parse.y parse.l
	flex parse.l
	bison -d parse.y 
	gcc -o test_parse parse.tab.c lex.yy.c -lfl

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ out.make.dep Makefile.bak *.unit $(TESTS)
	-@rm parse.c *.out parse.tab.? lex.yy.c

#-----------------------------------------------------------------------

.PHONY: install
install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)   
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/parameters.hpp $(INCLUDE)


#=======================================================================
#	INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include

# DO NOT DELETE
