#
# test         Run unit tests
# clean        Clean directory
# test_hdf5    Compile the unit tests
# dep          Generate dependencies
# 

include ../Make.include

OBJ        = test_hdf5.o hdf5.o
INC_LOCAL  = -I../data
OBJ_EXTERN = ../data/array.o

OBJS = $(OBJ) $(OBJ_EXTERN)

test_hdf5: $(OBJS)
	$(CPP) -o test_hdf5 $(OBJS) $(INC) $(INC_IO) -lhdf5

test: test_hdf5
	test_hdf5

.PHONY: clean
clean:
	-rm *.o test_hdf5 *~ *.hdf5 out.make.dep Makefile.bak

dep:
	-@(makedepend $(INC) -- -- *.?pp) >& out.make.dep
# DO NOT DELETE

hdf5.o: /usr/include/hdf5.h /usr/include/H5public.h /usr/include/H5pubconf.h
hdf5.o: /usr/include/features.h /usr/include/sys/cdefs.h
hdf5.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
hdf5.o: /usr/include/gnu/stubs-32.h /usr/include/sys/types.h
hdf5.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
hdf5.o: /usr/include/time.h /usr/include/endian.h /usr/include/bits/endian.h
hdf5.o: /usr/include/sys/select.h /usr/include/bits/select.h
hdf5.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
hdf5.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
hdf5.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
hdf5.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
hdf5.o: /usr/include/bits/posix2_lim.h /usr/include/stdint.h
hdf5.o: /usr/include/bits/wchar.h /usr/include/inttypes.h
hdf5.o: /usr/include/H5api_adpt.h /usr/include/H5Apublic.h
hdf5.o: /usr/include/H5Ipublic.h /usr/include/H5ACpublic.h
hdf5.o: /usr/include/H5Bpublic.h /usr/include/H5Dpublic.h
hdf5.o: /usr/include/H5Epublic.h /usr/include/stdio.h /usr/include/libio.h
hdf5.o: /usr/include/_G_config.h /usr/include/wchar.h
hdf5.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
hdf5.o: /usr/include/H5Fpublic.h /usr/include/H5FDpublic.h
hdf5.o: /usr/include/H5Gpublic.h /usr/include/H5Opublic.h
hdf5.o: /usr/include/H5HGpublic.h /usr/include/H5HLpublic.h
hdf5.o: /usr/include/H5MMpublic.h /usr/include/H5Ppublic.h
hdf5.o: /usr/include/H5Zpublic.h /usr/include/H5Rpublic.h
hdf5.o: /usr/include/H5Spublic.h /usr/include/H5Tpublic.h
hdf5.o: /usr/include/H5FDcore.h /usr/include/H5FDfamily.h
hdf5.o: /usr/include/H5FDgass.h /usr/include/string.h /usr/include/H5FDlog.h
hdf5.o: /usr/include/H5FDmpi.h /usr/include/H5FDmpio.h
hdf5.o: /usr/include/H5FDmpiposix.h /usr/include/H5FDmulti.h
hdf5.o: /usr/include/H5FDsec2.h /usr/include/H5FDsrb.h
hdf5.o: /usr/include/H5FDstdio.h /usr/include/H5FDstream.h
hdf5.o: ../include/error.hpp ../include/scalar.hpp hdf5.hpp
test_hdf5.o: /usr/include/stdio.h /usr/include/features.h
test_hdf5.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
test_hdf5.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
test_hdf5.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
test_hdf5.o: /usr/include/libio.h /usr/include/_G_config.h
test_hdf5.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
test_hdf5.o: /usr/include/bits/sys_errlist.h /usr/include/hdf5.h
test_hdf5.o: /usr/include/H5public.h /usr/include/H5pubconf.h
test_hdf5.o: /usr/include/sys/types.h /usr/include/time.h
test_hdf5.o: /usr/include/endian.h /usr/include/bits/endian.h
test_hdf5.o: /usr/include/sys/select.h /usr/include/bits/select.h
test_hdf5.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
test_hdf5.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
test_hdf5.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
test_hdf5.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
test_hdf5.o: /usr/include/bits/posix2_lim.h /usr/include/stdint.h
test_hdf5.o: /usr/include/bits/wchar.h /usr/include/inttypes.h
test_hdf5.o: /usr/include/H5api_adpt.h /usr/include/H5Apublic.h
test_hdf5.o: /usr/include/H5Ipublic.h /usr/include/H5ACpublic.h
test_hdf5.o: /usr/include/H5Bpublic.h /usr/include/H5Dpublic.h
test_hdf5.o: /usr/include/H5Epublic.h /usr/include/H5Fpublic.h
test_hdf5.o: /usr/include/H5FDpublic.h /usr/include/H5Gpublic.h
test_hdf5.o: /usr/include/H5Opublic.h /usr/include/H5HGpublic.h
test_hdf5.o: /usr/include/H5HLpublic.h /usr/include/H5MMpublic.h
test_hdf5.o: /usr/include/H5Ppublic.h /usr/include/H5Zpublic.h
test_hdf5.o: /usr/include/H5Rpublic.h /usr/include/H5Spublic.h
test_hdf5.o: /usr/include/H5Tpublic.h /usr/include/H5FDcore.h
test_hdf5.o: /usr/include/H5FDfamily.h /usr/include/H5FDgass.h
test_hdf5.o: /usr/include/string.h /usr/include/H5FDlog.h
test_hdf5.o: /usr/include/H5FDmpi.h /usr/include/H5FDmpio.h
test_hdf5.o: /usr/include/H5FDmpiposix.h /usr/include/H5FDmulti.h
test_hdf5.o: /usr/include/H5FDsec2.h /usr/include/H5FDsrb.h
test_hdf5.o: /usr/include/H5FDstdio.h /usr/include/H5FDstream.h
test_hdf5.o: ../include/unit.hpp ../include/scalar.hpp hdf5.hpp
