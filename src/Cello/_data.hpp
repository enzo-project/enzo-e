// See LICENSE_CELLO file for license and copyright information

/// @file     _data.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-17
/// @brief    Private include file for the \ref Data component

#ifndef _DATA_HPP
#define _DATA_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

extern void png_array (const char * filename,
		       float * array,
		       int gx,int gy,int gz,
		       int mx,int my,int mz,
		       const char * file, int line,
		       int axis = 2,
		       int px=0, int py=0,
		       double scale = 1.0
		       );
//----------------------------------------------------------------------
// Defines
//----------------------------------------------------------------------

#define PARTICLE_ALIGN 16

// integer limits on particle position within a Block:
//
//  -N    -N/2   0    N/2    N
//   *-----*-----*-----*-----*
//  LEFT   | IN BLOCK  | RIGHT
//        PMIN        PMAX
//   Full integer range is [-N,N)
//   Particles in block [-N/2,N/2)

#define PMIN_8 -64
#define PMAX_8  64

#define PMIN_16 -16384
#define PMAX_16  16384

#define PMIN_32 -1073741824
#define PMAX_32  1073741824

#define PMIN_64 -4611686018427387904
#define PMAX_64  4611686018427387904

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

class FieldData;
class FieldFace;

#include "data_Grouping.hpp"

#include "data_FieldDescr.hpp"
#include "data_FieldData.hpp"
#include "data_Field.hpp"
#include "data_FieldFace.hpp"

#include "data_ItIndex.hpp"
#include "data_ItIndexList.hpp"
#include "data_ItIndexRange.hpp"

#include "data_ParticleDescr.hpp"
#include "data_ParticleData.hpp"
#include "data_Particle.hpp"

#include "data_Data.hpp"

#include "data_DataMsg.hpp"

#ifdef DEBUG_FIELD

#   define TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,PLOT)	\
  {									\
    double sum_all=0.0;							\
    double sum_real=0.0;						\
    double sum2_all=0.0;						\
    double sum2_real=0.0;						\
    for (int iz=0; iz<mz; iz++) {					\
      for (int iy=0; iy<my; iy++) {					\
	for (int ix=0; ix<mx; ix++) {					\
	  int i = ix + mx*(iy + my*iz);					\
	  sum_all+=SCALE*(FIELD[i]);					\
	  sum2_all += SCALE*(FIELD[i]) * SCALE*(FIELD[i]);		\
	}								\
      }									\
    }									\
    for (int iz=gz; iz<mz-gz; iz++) {					\
      for (int iy=gy; iy<my-gy; iy++) {					\
	for (int ix=gx; ix<mx-gx; ix++) {				\
	  int i = ix + mx*(iy + my*iz);					\
	  if (PLOT) CkPrintf ("PRINT:%s %24.18f\n",NAME,SCALE*FIELD[i]);	\
	  sum_real+=SCALE*(FIELD[i]);					\
	  sum2_real+=SCALE*(FIELD[i])*SCALE*(FIELD[i]);			\
	}								\
      }									\
    }									\
    CkPrintf ("DEBUG_FIELD (%g) [%g] | (%g : %g %g %g) %s:%d %s\n",	\
	      sum_real,sum_all,					\
	      SCALE*(FIELD[gx+mx*(gy+my*gz)]),				\
	      SCALE*(FIELD[(gx+1)+mx*(gy+my*gz)]),			\
	      SCALE*(FIELD[gx+mx*((gy+1)+my*gz)]),			\
	      SCALE*(FIELD[gx+mx*(gy+my*(gz+1))]),			\
	      __FILE__,__LINE__,NAME);					\
    fflush(stdout);							\
    char filename[80];						\
    sprintf (filename,"renzo-p-%s.png",NAME);				\
    png_array (filename,(float*)(FIELD),gx,gy,gz,mx,my,mz,__FILE__,__LINE__,2,16,16,SCALE); \
    sprintf (filename,"enzo-p-%s.png",NAME);				\
    png_array (filename,(float*)(FIELD),0,0,0,mx,my,mz,__FILE__,__LINE__,2,16,16,SCALE); \
  }
#   define TRACE_FIELD_(NAME,FIELD,SCALE) TRACE_FIELD_GM(NAME,FIELD,SCALE,gx_,gy_,gz_,mx_,my_,mz_,false)
#   define TRACE_FIELD(NAME,FIELD,SCALE)  TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,false)
#   define PRINT_FIELD_(NAME,FIELD,SCALE) TRACE_FIELD_GM(NAME,FIELD,SCALE,gx_,gy_,gz_,mx_,my_,mz_,true)
#   define PRINT_FIELD(NAME,FIELD,SCALE)  TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,true)

#else

#   define TRACE_FIELD(NAME,FIELD,SCALE) /* ... */
#   define PRINT_FIELD(NAME,FIELD,SCALE) /* ... */
#   define TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz) /* ... */
#   define TRACE_FIELD_(NAME,FIELD,SCALE) /* ... */

#endif

#endif /* _DATA_HPP */

