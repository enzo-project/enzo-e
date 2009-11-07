TOP     = ../..
PACKAGE = disk
OBJS    = disk_hdf5.o disk_ifrit.o
TESTS   = test_disk_hdf5 test_disk_ifrit

#-----------------------------------------------------------------------

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

test_disk_hdf5: test_disk_hdf5.o 
	$(CXX) -o test_disk_hdf5 test_disk_hdf5.o $(LIB)
test_disk_ifrit: test_disk_ifrit.o 
	$(CXX) -o test_disk_ifrit test_disk_ifrit.o $(LIB)
#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/disk.hpp      $(INCLUDE)
	ln -sf `pwd`/disk_hdf5.hpp $(INCLUDE)
	ln -sf `pwd`/disk_ifrit.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
