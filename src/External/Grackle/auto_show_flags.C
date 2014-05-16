#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"CPP = \n");
   fprintf (fp,"CC  = /usr/bin/gcc\n");
   fprintf (fp,"CXX = /usr/bin/g++\n");
   fprintf (fp,"FC  = /usr/bin/gfortran\n");
   fprintf (fp,"F90 = /usr/bin/gfortran\n");
   fprintf (fp,"LD  = /usr/bin/g++\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API -fPIC -DLARGE_INTS  -DCONFIG_BFLOAT_8\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/usr/local/hdf5/1.8.2s/include     -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = -P -traditional\n");
   fprintf (fp,"CFLAGS   =  -O2\n");
   fprintf (fp,"CXXFLAGS =  -O2\n");
   fprintf (fp,"FFLAGS   = -fno-second-underscore -ffixed-line-length-132 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8  -O2\n");
   fprintf (fp,"F90FLAGS = -fno-second-underscore -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8  -O2\n");
   fprintf (fp,"LDFLAGS  =  -O2\n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/usr/local/hdf5/1.8.2s/lib -lhdf5 -lz  -lgfortran \n");
   fprintf (fp,"\n");
}
