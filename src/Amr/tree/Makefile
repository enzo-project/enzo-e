
TOP     = ../../..
PACKAGE = amr
OBJS    = amr_tree4.o  amr_node4.o \
          amr_tree16.o amr_node16.o \
          amr_tree2k.o amr_node2k.o
TESTS   = test_tree test_tree_k
DIRS    = 

#-----------------------------------------------------------------------

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

test_tree: test_tree.o $(OBJS)
	$(CXX) -o test_tree test_tree.o $(LIB)
test_tree_k: test_tree_k.o $(OBJS)
	$(CXX) -o test_tree_k test_tree_k.o $(LIB) #-lpng 

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/amr_node.hpp $(INCLUDE)
	ln -sf `pwd`/amr_tree4.hpp $(INCLUDE)
	ln -sf `pwd`/amr_node4.hpp $(INCLUDE)
	ln -sf `pwd`/amr_tree16.hpp $(INCLUDE)
	ln -sf `pwd`/amr_node16.hpp $(INCLUDE)
	ln -sf `pwd`/amr_tree2k.hpp $(INCLUDE)
	ln -sf `pwd`/amr_node2k.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
