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

#include "data_ScalarDescr.hpp"
#include "data_ScalarData.hpp"
#include "data_Scalar.hpp"

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
    int gx=3,gy=3,gz=3;							\
    double sum_all=0.0;							\
    double sum_real=0.0;						\
    double sum_abs=0.0;							\
    double sum_mean=0.0;						\
    double sum_var=0.0;							\
    double sum2_all=0.0;						\
    double sum2_real=0.0;						\
    for (int iz=0; iz<mz; iz++) {					\
      for (int iy=0; iy<my; iy++) {					\
	for (int ix=0; ix<mx; ix++) {					\
	  int i = ix + mx*(iy + my*iz);					\
	  double value=SCALE*FIELD[i];					\
	  sum_all+=value;						\
	  sum2_all += value * value;					\
	}								\
      }									\
    }									\
    for (int iz=gz; iz<mz-gz; iz++) {					\
      for (int iy=gy; iy<my-gy; iy++) {					\
	for (int ix=gx; ix<mx-gx; ix++) {				\
	  int i = ix + mx*(iy + my*iz);					\
	  double value=SCALE*FIELD[i];					\
	  sum_real+=value;						\
	  sum_abs+=std::abs(value);					\
	  sum2_real+=value*value;					\
	}								\
      }									\
    }									\
    double mean=sum_real/((mx-2*gx)*(my-2*gy)*(mz-2*gz));		\
    for (int iz=gz; iz<mz-gz; iz++) {					\
      for (int iy=gy; iy<my-gy; iy++) {					\
	for (int ix=gx; ix<mx-gx; ix++) {				\
	  int i = ix + mx*(iy + my*iz);					\
	  double value=SCALE*FIELD[i];					\
	  sum_mean +=std::abs(value-mean);				\
	  sum_var  +=(value-mean)*(value-mean);				\
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
    CkPrintf ("DEBUG_FIELD %s %20.18g %20.18g %20.18g\n",  NAME,sum_abs,sum_mean,sum_var); \
    fflush(stdout);							\
    char filename[80];						\
    sprintf (filename,"renzo-e-%s.png",NAME);				\
    png_array (filename,(float*)(FIELD),gx,gy,gz,mx,my,mz,__FILE__,__LINE__,2,16,16,SCALE); \
    sprintf (filename,"enzo-e-%s.png",NAME);				\
    png_array (filename,(float*)(FIELD),0,0,0,mx,my,mz,__FILE__,__LINE__,2,16,16,SCALE); \
  }
#   define TRACE_FIELD_(NAME,FIELD,SCALE) TRACE_FIELD_GM(NAME,FIELD,SCALE,gx_,gy_,gz_,mx_,my_,mz_,false)
#   define TRACE_FIELD(NAME,FIELD,SCALE)  TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,false)
#   define PRINT_FIELD_(NAME,FIELD,SCALE) TRACE_FIELD_GM(NAME,FIELD,SCALE,gx_,gy_,gz_,mx_,my_,mz_,true)
#   define PRINT_FIELD(NAME,FIELD,SCALE)  TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz,true)

#define TRACE_PARTICLE(ATTRIBUTE,particle,TYPE_NAME,ATTR_NAME)		\
  {									\
    int it = particle.type_index(TYPE_NAME);				\
    int ia = particle.attribute_index(it,ATTR_NAME);			\
    int nb = particle.num_batches(it);					\
    double sum_abs = 0.0;						\
    for (int ib=0; ib<nb; ib++) {					\
      int np = particle.num_particles(it,ib);				\
      enzo_float * array = (enzo_float *) particle.attribute_array(it,ia,ib); \
      for (int ip=0; ip<np; ip++) {					\
	sum_abs += std::abs(array[ip]);					\
      }									\
    }									\
    CkPrintf ("DEBUG_PARTICLE %s %20.18g\n",  ATTRIBUTE,sum_abs);		\
  }


#else

#   define TRACE_FIELD(NAME,FIELD,SCALE) /* ... */
#   define PRINT_FIELD(NAME,FIELD,SCALE) /* ... */
#   define TRACE_FIELD_GM(NAME,FIELD,SCALE,gx,gy,gz,mx,my,mz) /* ... */
#   define TRACE_FIELD_(NAME,FIELD,SCALE) /* ... */
#   define TRACE_PARTICLE(ATTRIBUTE,particle,TYPE_NAME,ATTR_NAME) /* ... */

#endif

#endif /* _DATA_HPP */

