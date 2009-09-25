
TOP     = ../..
PACKAGE = parameters
OBJS    = parameters.o
TESTS   = test_parameters test_parse

include $(TOP)/src/MAKE/Make.include

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

.PHONY: install-include
install-include:
	ln -sf `pwd`/parameters.hpp $(INCLUDE)
#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	-@(rm -f lex.yy.c parse.c parse.tab.? *.out)


