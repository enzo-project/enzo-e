SHELL = /bin/csh 

# Define MPCPP
include Make.machine

FLAGS = -g

INC = -I$(HOME)/include 
LIBS = -L$(HOME)/lib -lHYPRE_sstruct_ls -lHYPRE_sstruct_mv -lHYPRE_struct_ls -lHYPRE_struct_mv -lHYPRE_parcsr_ls -lHYPRE_parcsr_mv -lHYPRE_IJ_mv -lHYPRE_LSI -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_utilities
OBJS = \
       discret.o \
       grid.o \
       hierarchy.o \
       hypre-grav.o \
       hypre.o \
       level.o \
       mpi.o \
       point.o \
       problem.o \
       sphere.o

hypre-grav: dep $(OBJS)
hypre-grav: $(OBJS)
	$(MPCPP) -o hypre-grav $(OBJS) $(LIBS)


.cpp.o:
	$(MPCPP) -c -o $*.o $*.cpp $(FLAGS) $(INC)

.PHONY: dep
dep:
	-@(makedepend *.cpp *.hpp) >& /dev/null

.PHONY: clean
clean:
	rm -f *.o hypre-grav *~ \#*
doc:
	doxygen


# DO NOT DELETE

discret.o: /usr/include/assert.h /usr/include/features.h
discret.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
discret.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
discret.o: /usr/include/stdio.h /usr/include/bits/types.h
discret.o: /usr/include/bits/typesizes.h /usr/include/libio.h
discret.o: /usr/include/_G_config.h /usr/include/wchar.h
discret.o: /usr/include/bits/wchar.h /usr/include/gconv.h
discret.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
discret.o: /usr/include/math.h /usr/include/bits/huge_val.h
discret.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
discret.o: scalar.hpp discret.hpp mpi.hpp grid.hpp
grid.o: /usr/include/assert.h /usr/include/features.h
grid.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
grid.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
grid.o: /usr/include/stdio.h /usr/include/bits/types.h
grid.o: /usr/include/bits/typesizes.h /usr/include/libio.h
grid.o: /usr/include/_G_config.h /usr/include/wchar.h
grid.o: /usr/include/bits/wchar.h /usr/include/gconv.h
grid.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
grid.o: /usr/include/math.h /usr/include/bits/huge_val.h
grid.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h scalar.hpp
grid.o: discret.hpp mpi.hpp grid.hpp
hierarchy.o: /usr/include/assert.h /usr/include/features.h
hierarchy.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
hierarchy.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
hierarchy.o: /usr/include/stdio.h /usr/include/bits/types.h
hierarchy.o: /usr/include/bits/typesizes.h /usr/include/libio.h
hierarchy.o: /usr/include/_G_config.h /usr/include/wchar.h
hierarchy.o: /usr/include/bits/wchar.h /usr/include/gconv.h
hierarchy.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
hierarchy.o: scalar.hpp discret.hpp mpi.hpp grid.hpp level.hpp hierarchy.hpp
hypre.o: /usr/include/assert.h /usr/include/features.h
hypre.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
hypre.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h mpi.hpp
hypre.o: scalar.hpp point.hpp discret.hpp grid.hpp level.hpp hierarchy.hpp
hypre.o: sphere.hpp problem.hpp hypre.hpp
hypre-grav.o: /usr/include/assert.h /usr/include/features.h
hypre-grav.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
hypre-grav.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
hypre-grav.o: /usr/include/stdio.h /usr/include/bits/types.h
hypre-grav.o: /usr/include/bits/typesizes.h /usr/include/libio.h
hypre-grav.o: /usr/include/_G_config.h /usr/include/wchar.h
hypre-grav.o: /usr/include/bits/wchar.h /usr/include/gconv.h
hypre-grav.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
hypre-grav.o: /usr/include/stdlib.h /usr/include/sys/types.h
hypre-grav.o: /usr/include/time.h /usr/include/endian.h
hypre-grav.o: /usr/include/bits/endian.h /usr/include/sys/select.h
hypre-grav.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
hypre-grav.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
hypre-grav.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
hypre-grav.o: /usr/include/math.h /usr/include/bits/huge_val.h
hypre-grav.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
hypre-grav.o: scalar.hpp point.hpp sphere.hpp discret.hpp mpi.hpp grid.hpp
hypre-grav.o: level.hpp hierarchy.hpp problem.hpp hypre.hpp
level.o: /usr/include/assert.h /usr/include/features.h
level.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
level.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
level.o: /usr/include/stdio.h /usr/include/bits/types.h
level.o: /usr/include/bits/typesizes.h /usr/include/libio.h
level.o: /usr/include/_G_config.h /usr/include/wchar.h
level.o: /usr/include/bits/wchar.h /usr/include/gconv.h
level.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
level.o: scalar.hpp discret.hpp mpi.hpp grid.hpp level.hpp
mpi.o: /usr/include/stdio.h /usr/include/features.h /usr/include/sys/cdefs.h
mpi.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
mpi.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
mpi.o: /usr/include/bits/typesizes.h /usr/include/libio.h
mpi.o: /usr/include/_G_config.h /usr/include/wchar.h
mpi.o: /usr/include/bits/wchar.h /usr/include/gconv.h
mpi.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h mpi.hpp
point.o: /usr/include/assert.h /usr/include/features.h
point.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
point.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h scalar.hpp
point.o: point.hpp
problem.o: /usr/include/assert.h /usr/include/features.h
problem.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
problem.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
problem.o: /usr/include/stdio.h /usr/include/bits/types.h
problem.o: /usr/include/bits/typesizes.h /usr/include/libio.h
problem.o: /usr/include/_G_config.h /usr/include/wchar.h
problem.o: /usr/include/bits/wchar.h /usr/include/gconv.h
problem.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
problem.o: scalar.hpp point.hpp discret.hpp mpi.hpp grid.hpp level.hpp
problem.o: hierarchy.hpp sphere.hpp problem.hpp
sphere.o: /usr/include/assert.h /usr/include/features.h
sphere.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
sphere.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
sphere.o: /usr/include/stdio.h /usr/include/bits/types.h
sphere.o: /usr/include/bits/typesizes.h /usr/include/libio.h
sphere.o: /usr/include/_G_config.h /usr/include/wchar.h
sphere.o: /usr/include/bits/wchar.h /usr/include/gconv.h
sphere.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
sphere.o: scalar.hpp sphere.hpp
test.o: /usr/include/stdio.h /usr/include/features.h /usr/include/sys/cdefs.h
test.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
test.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
test.o: /usr/include/bits/typesizes.h /usr/include/libio.h
test.o: /usr/include/_G_config.h /usr/include/wchar.h
test.o: /usr/include/bits/wchar.h /usr/include/gconv.h
test.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
