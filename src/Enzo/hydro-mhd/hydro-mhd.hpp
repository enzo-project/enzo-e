// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/hydro-mhd/hydro-mhd.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for hydro-mhd subcomponent within the \ref Enzo layer

#ifndef ENZO_HYDROMHD_HYDROMHD_HPP
#define ENZO_HYDROMHD_HYDROMHD_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <array>
#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Method
#include "Cello/view.hpp" // CelloView

#include "Enzo/enzo.hpp" // enzo_float, EnzoConfig, EFltArrayMap

// in the future, the following should be removed from this header file
// (there's no NEED for it to be a transitive dependency for anything that
//  depends on the hydro-mhd dependency)
// - this currently to be included after the headers for EnzoEFltArrayMap,
//   EnzoCenteredFieldRegistry, & EnzoEOSVariant
// - TODO: make the Riemann Header Self-contained
#include "Enzo/hydro-mhd/riemann/EnzoRiemann.hpp"

//----------------------------------------------------------------------
// Assorted Public Functions
//----------------------------------------------------------------------

// TODO: move these into their own private header file
extern "C" void FORTRAN_NAME(calc_dt)
  (int *rank, int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   enzo_float *dx, enzo_float *dy, enzo_float *dz, 
   enzo_float *gamma, int *ipfree, enzo_float *aye,
   enzo_float *d, enzo_float *p, enzo_float *u, enzo_float *v, enzo_float *w,
   enzo_float *dt);

extern "C" void FORTRAN_NAME(calc_dt_ppml)
  (int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   enzo_float *dx, enzo_float *dy, enzo_float *dz,
   enzo_float *dn, enzo_float *vx, enzo_float *vy, enzo_float *vz, 
   enzo_float *bx, enzo_float *by, enzo_float *bz, 
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

//----------------------------------------------------------------------

#include "Enzo/hydro-mhd/toolkit/EnzoIntegrationQuanUpdate.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoLazyPassiveScalarFieldList.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoSourceGravity.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoSourceInternalEnergy.hpp"

// [order dependencies:]
#include "Enzo/hydro-mhd/toolkit/EnzoReconstructor.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoReconstructorNN.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoReconstructorPLM.hpp"

// [order dependencies:]
#include "Enzo/hydro-mhd/toolkit/EnzoBfieldMethod.hpp"
#include "Enzo/hydro-mhd/toolkit/EnzoBfieldMethodCT.hpp"

#include "Enzo/hydro-mhd/EnzoMethodMHDVlct.hpp"
#include "Enzo/hydro-mhd/EnzoMethodPpm.hpp"
#include "Enzo/hydro-mhd/EnzoMethodPpml.hpp"

#endif /* ENZO_HYDROMHD_HYDROMHD_HPP */
