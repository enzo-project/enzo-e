// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/hydro-mhd/ppm_fortran.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Function prototypes for ppm subroutines defined in Fortran

#ifndef ENZO_HYDROMHD_PPMFORTRAN_PPMFORTRAN_HPP
#define ENZO_HYDROMHD_PPMFORTRAN_PPMFORTRAN_HPP

#include "Enzo/enzo.hpp" // enzo_float, FORTRAN_NAME

extern "C" void FORTRAN_NAME(calc_dt)
  (int *rank, int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   enzo_float *dx, enzo_float *dy, enzo_float *dz,
   enzo_float *gamma, int *ipfree, enzo_float *aye,
   enzo_float *d, enzo_float *p, enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *dt);

extern "C" void FORTRAN_NAME(ppm_de)
  (enzo_float *d, enzo_float *E, enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *ge,
   int *grav, enzo_float *gr_ax, enzo_float *gr_ay, enzo_float *gr_az,
   enzo_float *gamma, enzo_float *dt, int *cycle_number,
   enzo_float dx[], enzo_float dy[], enzo_float dz[],
   int *rank, int *in, int *jn, int *kn,
   int is[], int ie[],
   int *flatten, int *ipresfree,
   int * iconsrec, int *iposrec,
   int *diff, int *steepen, int *idual,
   enzo_float *eta1, enzo_float *eta2,
   int *num_subgrids, int leftface[], int rightface[],
   int istart[], int iend[], int jstart[], int jend[],
   enzo_float *standard, int dindex[], int Eindex[],
   int uindex[], int vindex[], int windex[],
   int geindex[], enzo_float *temp,
   int *ncolor, enzo_float *colorpt, int *coloff,
   int colindex[], enzo_float *pressure_floor, enzo_float *density_floor, int *error,
   int *ie_error_x,
   int *ie_error_y,
   int *ie_error_z,
   int *num_ie_error
   );

#endif /* ENZO_HYDROMHD_PPMFORTRAN_PPMFORTRAN_HPP */
