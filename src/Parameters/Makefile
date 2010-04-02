
TOP     = ../..
PACKAGE = parameters
OBJS_PARSE = parse.tab.o lex.yy.o
OBJS    = parameters.o param.o $(OBJS_PARSE)
TESTS   = test_parameters
EXES    = test_parse

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: run-exes

run-exes: $(EXES)
	-./test_parse < implosion.in    > implosion.out
	-./test_parse < shock-pool.in   > shock-pool.out
	-./test_parse < test_groups.in  > test_groups.out
	-./test_parse < test_integer.in > test_integer.out
	-./test_parse < test.in         > test.out
	-./test_parse < test_list.in    > test_list.out
	-./test_parse < test_logical.in > test_logical.out
	-./test_parse < test_scalar.in  > test_scalar.out
	-grep -i -l error *.out

#-----------------------------------------------------------------------

test_parameters: test_parameters.o $(OBJS)
	$(CXX) -o test_parameters test_parameters.o $(LIB)

test_parse: test_parse.o $(OBJS_PARSE)
	$(CXX) -o test_parse test_parse.o $(LIB)

parse.tab.o:parse.tab.c
	$(CC) -o parse.tab.o parse.tab.c -c $(OPT)
lex.yy.c:parse.l
	flex parse.l
parse.tab.h: parse.y
	bison -d parse.y 
parse.tab.c: parse.y
	bison -d parse.y 
lex.yy.o:lex.yy.c
	$(CC) -o lex.yy.o lex.yy.c -c $(OPT)

.PHONY: install-include
install-include:
	ln -sf `pwd`/parameters.hpp $(INCLUDE)
	ln -sf `pwd`/param.hpp      $(INCLUDE)
	ln -sf `pwd`/parse.h        $(INCLUDE)
	ln -sf `pwd`/parse.tab.h    $(INCLUDE)
	ln -sf `pwd`/parameters_parameters.hpp $(INCLUDE)
#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	-@(rm -f lex.yy.c parse.c parse.tab.? *.out)


