// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/hydro-mhd/ppml_fortran.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Function prototypes for ppml subroutines defined in Fortran

#ifndef ENZO_HYDROMHD_PPMLFORTRAN_PPMLFORTRAN_HPP
#define ENZO_HYDROMHD_PPMLFORTRAN_PPMLFORTRAN_HPP

#include "Enzo/enzo.hpp" // enzo_float, FORTRAN_NAME

extern "C" void FORTRAN_NAME(calc_dt_ppml)
  (int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   enzo_float *dx, enzo_float *dy, enzo_float *dz,
   enzo_float *dn, enzo_float *vx, enzo_float *vy, enzo_float *vz,
   enzo_float *bx, enzo_float *by, enzo_float *bz,
   enzo_float *dt);

extern "C" void FORTRAN_NAME(ppml)
  (enzo_float *dn,   enzo_float *vx,   enzo_float *vy,   enzo_float *vz,
   enzo_float *bx,   enzo_float *by,   enzo_float *bz,
   enzo_float *dnrx, enzo_float *vxrx, enzo_float *vyrx, enzo_float *vzrx,
   enzo_float *bxrx, enzo_float *byrx, enzo_float *bzrx,
   enzo_float *dnry, enzo_float *vxry, enzo_float *vyry, enzo_float *vzry,
   enzo_float *bxry, enzo_float *byry, enzo_float *bzry,
   enzo_float *dnrz, enzo_float *vxrz, enzo_float *vyrz, enzo_float *vzrz,
   enzo_float *bxrz, enzo_float *byrz, enzo_float *bzrz,
   enzo_float *dt, enzo_float dx[], enzo_float dy[], enzo_float dz[],
   int *in, int *jn, int *kn,
   int is[], int ie[],
   int *num_subgrids, int leftface[], int rightface[],
   int istart[], int iend[], int jstart[], int jend[],
   enzo_float *standard, int dnindex[],
   int vxindex[], int vyindex[], int vzindex[],
   int bxindex[], int byindex[], int bzindex[],
   enzo_float *f1,enzo_float *f2,enzo_float *f3,enzo_float *f4,
   enzo_float *f5,enzo_float *f6,enzo_float *f7,
   enzo_float *g1,enzo_float *g2,enzo_float *g3,enzo_float *g4,
   enzo_float *g5,enzo_float *g6,enzo_float *g7,
   enzo_float *h1,enzo_float *h2,enzo_float *h3,enzo_float *h4,
   enzo_float *h5,enzo_float *h6,enzo_float *h7,
   enzo_float *ex,enzo_float *ey,enzo_float *ez,
   enzo_float *qu1,enzo_float *qu2,enzo_float *qu3,enzo_float *qu4,
   enzo_float *qu5,enzo_float *qu6,enzo_float *qu7);

#endif /* ENZO_HYDROMHD_PPMFORTRAN_PPMFORTRAN_HPP */
