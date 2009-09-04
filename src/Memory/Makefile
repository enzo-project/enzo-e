#
# (x) clean:        Clean directory
# (^) dep:          Update dependencies
# (X) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

PACKAGE   = memory
OBJ       = memory.o test_memory.o
TEST_EXES = \
        test_memory

#=======================================================================
OBJ_EXTERN = 

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: $(OBJS) $(TEST_EXES)

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

.PHONY: test
test: $(TEST_EXES) run-tests

#-----------------------------------------------------------------------

.PHONY: run-tests
run-tests: $(TEST_EXES)
	-./test_memory

#-----------------------------------------------------------------------

test_memory: test_memory.o 
	$(CXX) -o test_memory test_memory.o $(LIB) $(INC)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ out.make.dep Makefile.bak $(TEST_EXE)

#-----------------------------------------------------------------------

.PHONY: install
install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)   
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/memory.hpp              $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include

# DO NOT DELETE
