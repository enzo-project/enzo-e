# (X) install:         Install lib, include, and bin files
# (X) install-install: Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

TEST_EXES = \
	test_array \
	test_disk_hdf5 \
	test_error \
	test_parameters \
	test_performance_timer

TEST_OBJS = $(patsubst %,%.o,$(TEST_EXES))

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: test

#-----------------------------------------------------------------------

.PHONY: test
test: $(TEST_EXES)

#-----------------------------------------------------------------------

test_array: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_array $(TEST_OBJS) $(OBJS) $(INC)
test_performance_timer: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_performance_timer $(TEST_OBJS) $(OBJS) $(OPT) $(INC)
test_parameters: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_parameters $(TEST_OBJS) $(OBJS) $(INC)
test_error: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_error $(TEST_OBJS) $(OBJS) $(OPT) $(INC)
test_disk_hdf5: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_disk_hdf5 $(TEST_OBJS) $(OBJS) $(INC) $(LIB_HDF5)


test_test:

#-----------------------------------------------------------------------

libtest.a:

#-----------------------------------------------------------------------

.PHONY: clean
clean:

#-----------------------------------------------------------------------

.PHONY: install
install: $(TEST_EXES)
	cp ./test_performance_timer $(TOP)/bin
	cp ./test_array             $(TOP)/bin
	cp ./test_parameters        $(TOP)/bin
	cp ./test_error             $(TOP)/bin
	cp ./test_disk_hdf5         $(TOP)/bin


#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/test.hpp         $(INCLUDE)
	ln -sf `pwd`/test_unit_.hpp   $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include
