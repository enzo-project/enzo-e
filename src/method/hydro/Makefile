#
# test         Run unit tests
# clean        Clean directory
# test_array   Compile the unit tests
# dep          Generate dependencies
# 

include ../../Make.include

OBJS = 

test_hydro: $(OBJS)
	$(CPP) -o test_hydro $(OBJS) $(INC)

test: test_hydro
	test_hydro

.PHONY: clean
clean:
	-rm *.o test_hydro *~ out.make.dep Makefile.bak

dep:
	-@(makedepend $(INC) -- -- *.?pp) >& out.make.dep
# DO NOT DELETE

