TOP     = ../..
PACKAGE = disk
OBJS    = disk_hdf5.o
TESTS   = test_disk_hdf5

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

test_disk_hdf5: test_disk_hdf5.o 
	$(CXX) -o test_disk_hdf5 test_disk_hdf5.o $(LIB)
#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/disk.hpp      $(INCLUDE)
	ln -sf `pwd`/disk_hdf5.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
