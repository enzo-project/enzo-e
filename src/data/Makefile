#
# test         Run unit tests
# clean        Clean directory
# test_array   Compile the unit tests
# dep          Generate dependencies
# 

include ../Make.include

OBJS = test_array.o array.o

test_array: $(OBJS)
	$(CPP) -o test_array $(OBJS) $(INC)

test: test_array
	test_array

.PHONY: clean
clean:
	-rm *.o test_array *~ out.make.dep Makefile.bak

dep:
	-@(makedepend $(INC) -- -- *.?pp) >& out.make.dep
# DO NOT DELETE

array.o: /usr/include/stdio.h /usr/include/features.h
array.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
array.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
array.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
array.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
array.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
array.o: /usr/include/assert.h ../include/scalar.hpp ../include/error.hpp
array.o: array.hpp
test_array.o: /usr/include/stdio.h /usr/include/features.h
test_array.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
test_array.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
test_array.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
test_array.o: /usr/include/libio.h /usr/include/_G_config.h
test_array.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
test_array.o: /usr/include/bits/sys_errlist.h ../include/scalar.hpp array.hpp
test_array.o: ../include/unit.hpp /usr/include/string.h
