#
# test         Run unit tests
# clean        Clean directory
# test_array   Compile the unit tests
# dep          Generate dependencies
# 

TOP = ..

include $(TOP)/Make.include

TEST_OBJS = test_array.o 
OBJS = array.o

test_array: $(TEST_OBJS) libarray.a
	$(CPP) -o test_array $(TEST_OBJS) $(INC) -L. -larray
libarray.a: $(OBJS)
	ar rv $*.a $(OBJS)

test: test_array
	test_array

.PHONY: clean
clean:
	-rm *.o test_array *~ out.make.dep Makefile.bak libarray.a Array.unit

.PHONY: dep
dep:
	-@(makedepend -fDEPEND $(INC) -- -- *.?pp) >& out.make.dep

include DEPEND

