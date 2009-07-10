#
# (x) test:            Run unit tests
# (x) test_parameters: Compile the unit tests
# (x) clean:           Clean directory
# (^) dep:             Update dependencies
# (X) install:         Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

PACKAGE   = parameters
TEST_OBJS = test_parameters.o 
OBJ       = parameters.o

#=======================================================================
OBJ_EXTERN = 

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================
# TARGETS
#=======================================================================

.PHONY: test
test: test_parameters
	./test_parameters

#-----------------------------------------------------------------------

test_parameters: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_parameters $(TEST_OBJS) $(OBJS) $(INC)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o test_parameters *~ out.make.dep Makefile.bak *.unit

#-----------------------------------------------------------------------

.PHONY: install
install: install-include test_parameters $(OBJS)
	cp ./test_parameters $(TOP)/bin
	ar rcu $(TOP)/lib/libcello.a $(OBJS)   
	ar rcu $(TOP)/lib/libcello-$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/parameters.hpp               $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include

# DO NOT DELETE
