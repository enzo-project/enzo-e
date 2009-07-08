#
# (x) test:         Run unit tests
# (x) test_array:   Compile the unit tests
# (x) clean:        Clean directory
# (^) dep:          Update dependencies
# (X) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

TEST_OBJS = test_array.o 
OBJ       = array.o

#=======================================================================
INC_LOCAL  = -I$(SRC)/error
OBJ_EXTERN = 

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================
# TARGETS
#=======================================================================

.PHONY: test
test: test_array
	./test_array

#-----------------------------------------------------------------------

test_array: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_array $(TEST_OBJS) $(OBJS) $(INC)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o test_array *~ out.make.dep Makefile.bak *.unit

#-----------------------------------------------------------------------

.PHONY: install
install: install-include test_array $(OBJS)
	cp ./test_array $(TOP)/bin
	ar rcu $(TOP)/lib/libcello.a $(OBJS)   

#-----------------------------------------------------------------------

INSTALL = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/array.hpp         $(INCLUDE)
	ln -sf `pwd`/array_array_.hpp  $(INCLUDE)
	ln -sf `pwd`/array_scalar_.hpp $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include

# DO NOT DELETE
