
TOP = ../../..

TEST_EXES = test_tree

INC_LOCAL  = -I$(TOP)/include
LIB_LOCAL  = -L$(TOP)/lib -lcello


OBJS = tree4.o node4.o

FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3

FFLAGS   = $(FFLAGS_LOCAL)
CXXFLAGS = $(CXXFLAGS_LOCAL) $(INC_LOCAL)
LIB      = $(LIB_LOCAL)  -lhdf5

all: $(TEST_EXES)
test_tree: $(OBJS) test_tree.o
	$(CXX) -o test_tree test_tree.o $(OBJS) $(LIB)

.PHONY: install
install: $(TEST_EXES)
	cp $(TEST_EXES) $(TOP)/bin

obj: $(OBJS)

test: $(TEST_EXES)
#	./test_implosion # slow
#	./test_enzo      # slow

.PHONY: clean
clean:
	-@rm *.o $(TEST_EXES) *~

include $(TOP)/src/MAKE/Make.include

include DEPEND
