// See LICENSE_CELLO file for license and copyright information

/// @file     field_RestrictLinear.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-10
/// @brief    Implentation of default linear restriction

#include "problem.hpp"

//----------------------------------------------------------------------

RestrictLinear::RestrictLinear() throw()
  : Restrict ()
{
  TRACE("RestrictLinear::RestrictLinear");
}

//----------------------------------------------------------------------

void RestrictLinear::apply 
( precision_type precision,
  void *       values_c, int nd3_c[3], int n3_c[3],
  const void * values_f, int nd3_f[3], int n3_f[3])
{
  switch (precision) {

  case precision_single:

    apply_( (float *) values_c,       nd3_c, n3_c,
	    (const float *) values_f, nd3_f, n3_f);

    break;

  case precision_double:

    apply_( (double *) values_c,       nd3_c, n3_c,
	    (const double *) values_f, nd3_f, n3_f);

    break;

  default:

    ERROR1 ("RestrictLinear::apply()",
	    "Unknown precision %d",
	    precision);
    break;
  }
}

//----------------------------------------------------------------------

template<class T>
void RestrictLinear::apply_
( T *       values_c, int nd3_c[3], int n3_c[3],
  const T * values_f, int nd3_f[3], int n3_f[3])
{
  const int rank = (nd3_f[1] == 1) ? 1 : ((nd3_f[2] == 1) ? 2 : 3);

  const int dx = 1;
  const int dy = nd3_f[0];
  const int dz = nd3_f[0]*nd3_f[1];

  if (rank == 1) {
    for (int i_c=0; i_c<n3_c[0]; i_c++) {
      int i_f = i_c*2;
      values_c[i_c] = 
	values_f[i_f     ] + 
	values_f[i_f + dx];
    }

  } else if (rank == 2) {
    for (int ix_c=0; ix_c<n3_c[0]; ix_c++) {
      int ix_f = ix_c*2;
      for (int iy_c=0; iy_c<n3_c[1]; iy_c++) {
	int iy_f = iy_c*2;
	int i_c = ix_c + nd3_c[0]*iy_c;
	int i_f = ix_f + nd3_f[0]*iy_f;
	values_c[ix_c] = 
	  values_f[i_f          ] + 
	  values_f[i_f + dx     ] +
	  values_f[i_f +      dy] + 
	  values_f[i_f + dx + dy];

      }
    }

  } else if (rank == 3) {
    for (int ix_c=0; ix_c<n3_c[0]; ix_c++) {
      int ix_f = ix_c*2;
      for (int iy_c=0; iy_c<n3_c[1]; iy_c++) {
	int iy_f = iy_c*2;
	for (int iz_c=0; iz_c<n3_c[2]; iz_c++) {
	  int iz_f = iz_c*2;
	  int i_c = ix_c + nd3_c[0]*iy_c;
	  int i_f = ix_f + nd3_f[0]*iy_f;
	  values_c[ix_c] = 
	    values_f[i_f               ] + 
	    values_f[i_f + dx          ] +
	    values_f[i_f +      dy     ] + 
	    values_f[i_f + dx + dy     ] +
	    values_f[i_f           + dz] + 
	    values_f[i_f + dx      + dz] +
	    values_f[i_f +      dy + dz] + 
	    values_f[i_f + dx + dy + dz];
	}
      }
    }
  }
}

//======================================================================

