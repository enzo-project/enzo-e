#
# ( ) test:         Run unit tests
# ( ) test_error:   Compile the unit tests
# ( ) liberror.a:
# ( ) clean:        Clean directory
# ( ) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

TEST_OBJS = test_error.o
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

all: test

#-----------------------------------------------------------------------

.PHONY: test

test: test_error
	./test_error

#-----------------------------------------------------------------------

test_error: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_error $(TEST_OBJS) $(OBJS) $(OPT) $(INC)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ test_error lib*.a *.unit

#-----------------------------------------------------------------------

.PHONY: install

install: install-include test_error $(OBJS)
	cp ./test_error $(TOP)/bin
	ar rcu $(TOP)/lib/libcello.a $(OBJS)

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/error.hpp  $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include




