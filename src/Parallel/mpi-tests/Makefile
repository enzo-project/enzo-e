TOP     = ../../..
PACKAGE = parallel
OBJS    = mpi-distance.o mpi-array.o
TESTS   = mpi-distance mpi-array

include $(TOP)/src/MAKE/Make.include

#-----------------------------------------------------------------------

mpi-distance: mpi-distance.cpp
	$(CXX) -o mpi-distance mpi-distance.o $(LIB) $(LIB_HDF5)
mpi-array: mpi-array.cpp
	$(CXX) -o mpi-array mpi-array.o $(LIB) $(LIB_HDF5)

.PHONY: install-include
install-include:

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
