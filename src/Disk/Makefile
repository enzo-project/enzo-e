#
# (x) test:         Run unit tests
# (x) test_array:   Compile the unit tests
# (x) clean:        Clean directory
# (^) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

PACKAGE   = disk
TEST_OBJS = test_hdf5.o 
OBJ       = hdf5.o

#=======================================================================
OBJ_EXTERN = $(patsubst %,$(SRC)/%,$(ARRAY_OBJS))

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================

#=======================================================================
# TARGETS
#=======================================================================

#-----------------------------------------------------------------------

.PHONY: test
test: test_hdf5
	./test_hdf5

#-----------------------------------------------------------------------

test_hdf5: $(TEST_OBJS) $(OBJS)
	echo test $(TEST_OBJS) $(OBJS) $(ARRAY_OBJ)
	$(CXX) -o test_hdf5 $(TEST_OBJS) $(OBJS) $(INC) $(LIB_HDF5)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o test_*hdf5 *~ out.make.dep Makefile.bak lib*.a *.unit *.hdf5

#-----------------------------------------------------------------------

.PHONY: install
install: install-include test_hdf5 $(OBJS)
	cp ./test_hdf5 $(TOP)/bin
	ar rcu $(TOP)/lib/libcello.a $(OBJS)
	ar rcu $(TOP)/lib/libcello-$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/disk.hpp       $(INCLUDE)
	ln -sf `pwd`/disk_hdf5.hpp $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include
include $(SRC)/MAKE/Make.objects



