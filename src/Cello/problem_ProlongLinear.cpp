// See LICENSE_CELLO file for license and copyright information

/// @file     problem_ProlongLinear.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of default linear prolongation

#include "problem.hpp"

//----------------------------------------------------------------------

ProlongLinear::ProlongLinear() throw()
  : Prolong ()
{
  TRACE("ProlongLinear::ProlongLinear");
}

//----------------------------------------------------------------------

int ProlongLinear::apply 
( precision_type precision,
  void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
  const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3])
{
  TRACE6("ProlongLinear fine   %d:%d %d:%d %d:%d",
         im3_f[0],n3_f[0]+im3_f[0],
         im3_f[1],n3_f[1]+im3_f[1],
         im3_f[2],n3_f[2]+im3_f[2]);

  TRACE6("ProlongLinear coarse %d:%d %d:%d %d:%d",
         im3_c[0],n3_c[0]+im3_c[0],
         im3_c[1],n3_c[1]+im3_c[1],
         im3_c[2],n3_c[2]+im3_c[2]);

  switch (precision)  {

  case precision_single:

    return apply_(       (float *) values_f, nd3_f, im3_f, n3_f,
			 (const float *) values_c, nd3_c, im3_c, n3_c);

    break;

  case precision_double:

    return apply_(       (double *) values_f, nd3_f, im3_f, n3_f,
			 (const double *) values_c, nd3_c, im3_c, n3_c);

    break;

  default:

    ERROR1 ("ProlongLinear::apply()",
            "Unknown precision %d",
            precision);

    return 0;
  }
}

//----------------------------------------------------------------------

template <class T>
int ProlongLinear::apply_
(       T * values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
	const T * values_c, int nd3_c[3], int im3_c[3], int n3_c[3])
{

  const int dx_c = 1;
  const int dy_c = nd3_c[0];
  const int dz_c = nd3_c[0]*nd3_c[1];

  int rank = (nd3_f[2] > 1) ? 3 : ( (nd3_f[1] > 1) ? 2 : 1 );

  for (int i=0; i<rank; i++) {
    const char * xyz = "xyz";
    ASSERT3 ("ProlongLinear::apply_",
             "fine array %c-axis %d must be 2 x coarse axis %d",
             xyz[i],n3_c[i],n3_f[i],
             n3_f[i]==n3_c[i]*2);

    ASSERT2 ("ProlongLinear::apply_",
             "fine grid %c-axis %d must be divisible by 2",
             xyz[i],n3_f[i],n3_f[i] % 2 == 0);
    ASSERT2 ("ProlongLinear::apply_",
             "fine grid %c-axis %d must be at least 4",
             xyz[i],n3_f[i],n3_f[i] >= 4);
  }

  if (n3_f[1]==1) {

    const int ix0_f = im3_f[0];
    const int ix0_c = im3_c[0];
    const int nx_f = n3_f[0];

    for (int ix_f = 0; ix_f<nx_f; ix_f++) {

      int ix_c = ((ix_f+1) >> 1) - 1;

      // shift coarse cells for fine grid cells on faces
      ix_c += (ix_f==0) ? +1 : (ix_f < nx_f-1) ? 0 : -1;

      // weighting factor for even / odd odd / even
      // take into account extrapolating at edges
      //
      // 1D: (5/4 -1/4) (3/4 1/4) (1/4 3/4) ... (3/4 1/4) (1/4 3/4) (-1/4 5/4)
      //     so wy is (5/4 -1/4) or (1/4 3/4) for even,
      //             (-1/4 5/4)  or (3/4 1/4) for odd
      
      int wx[2] = { (0 < ix_f && ix_f < nx_f-1) ? 1 : 5,
		    (0 < ix_f && ix_f < nx_f-1) ? 3 : -1 };
    
      T wx0 = 0.25*wx[ ix_f&1];
      T wx1 = 0.25*wx[~ix_f&1];

      int i_c = (ix0_c+ix_c) ;
      int i_f = (ix0_f+ix_f) ;
	  
      values_f[i_f] = wx0*values_c[i_c]
	+             wx1*values_c[i_c + dx_c ];
    }

    return (sizeof(T) * n3_c[0]);


  } else if (n3_f[2] == 1) {

    const int ix0_f = im3_f[0];
    const int iy0_f = im3_f[1];
    const int ix0_c = im3_c[0];
    const int iy0_c = im3_c[1];
    const int nx_f = n3_f[0];
    const int ny_f = n3_f[1];
    const int mx_c = nd3_c[0];
    const int mx_f = nd3_f[0];
    
    for (int iy_f = 0; iy_f<ny_f; iy_f++) {

      int iy_c = ((iy_f+1) >> 1) - 1;
      iy_c += (iy_f==0) ? +1 : (iy_f < ny_f-1) ? 0 : -1;

      int wy[2] = { (0 < iy_f && iy_f < ny_f-1) ? 1 : 5,
		    (0 < iy_f && iy_f < ny_f-1) ? 3 : -1 };
    
      T wy0 = 0.25*wy[ iy_f&1];
      T wy1 = 0.25*wy[~iy_f&1];

      for (int ix_f = 0; ix_f<nx_f; ix_f++) {

	int ix_c = ((ix_f+1) >> 1) - 1;
	ix_c += (ix_f==0) ? +1 : (ix_f < nx_f-1) ? 0 : -1;

	int wx[2] = { (0 < ix_f && ix_f < nx_f-1) ? 1 : 5,
		      (0 < ix_f && ix_f < nx_f-1) ? 3 : -1 };

	T wx0 = 0.25*wx[ ix_f&1];
	T wx1 = 0.25*wx[~ix_f&1];

	int i_c = (ix0_c+ix_c) + mx_c * ( (iy0_c+iy_c) );
	int i_f = (ix0_f+ix_f) + mx_f * ( (iy0_f+iy_f) );

	values_f[i_f] = wx0*wy0*values_c[i_c]
	  +             wx1*wy0*values_c[i_c + dx_c ]
	  +             wx0*wy1*values_c[i_c        + dy_c ]
	  +             wx1*wy1*values_c[i_c + dx_c + dy_c ];
      }
    }

    return (sizeof(T) * n3_c[0]*n3_c[1]);

  } else {

    const int ix0_f = im3_f[0];
    const int iy0_f = im3_f[1];
    const int iz0_f = im3_f[2];
    const int ix0_c = im3_c[0];
    const int iy0_c = im3_c[1];
    const int iz0_c = im3_c[2];
    const int nx_f = n3_f[0];
    const int ny_f = n3_f[1];
    const int nz_f = n3_f[2];
    const int mx_c = nd3_c[0];
    const int mx_f = nd3_f[0];
    const int my_c = nd3_c[1];
    const int my_f = nd3_f[1];
    
    for (int iz_f = 0; iz_f<nz_f; iz_f++) {

      int iz_c = ((iz_f+1) >> 1) - 1;
      iz_c += (iz_f==0) ? +1 : (iz_f < nz_f-1) ? 0 : -1;

      int wz[2] = { (0 < iz_f && iz_f < nz_f-1) ? 1 : 5,
		    (0 < iz_f && iz_f < nz_f-1) ? 3 : -1 };
    
      T wz0 = 0.25*wz[ iz_f&1];
      T wz1 = 0.25*wz[~iz_f&1];

      for (int iy_f = 0; iy_f<ny_f; iy_f++) {

	int iy_c = ((iy_f+1) >> 1) - 1;
	iy_c += (iy_f==0) ? +1 : (iy_f < ny_f-1) ? 0 : -1;

	int wy[2] = { (0 < iy_f && iy_f < ny_f-1) ? 1 : 5,
		      (0 < iy_f && iy_f < ny_f-1) ? 3 : -1 };

	T wy0 = 0.25*wy[ iy_f&1];
	T wy1 = 0.25*wy[~iy_f&1];

	for (int ix_f = 0; ix_f<nx_f; ix_f++) {

	  int ix_c = ((ix_f+1) >> 1) - 1;
	  ix_c += (ix_f==0) ? +1 : (ix_f < nx_f-1) ? 0 : -1;

	  int wx[2] = { (0 < ix_f && ix_f < nx_f-1) ? 1 : 5,
			(0 < ix_f && ix_f < nx_f-1) ? 3 : -1 };

	  T wx0 = 0.25*wx[ ix_f&1];
	  T wx1 = 0.25*wx[~ix_f&1];

	  int i_c = (ix0_c+ix_c) + mx_c*( (iy0_c+iy_c) + my_c*(iz0_c+iz_c) );
	  int i_f = (ix0_f+ix_f) + mx_f*( (iy0_f+iy_f) + my_f*(iz0_f+iz_f) );
	  
	  values_f[i_f] = wx0*wy0*wz0*values_c[i_c]
	    +             wx1*wy0*wz0*values_c[i_c + dx_c ]
	    +             wx0*wy1*wz0*values_c[i_c        + dy_c ]
	    +             wx1*wy1*wz0*values_c[i_c + dx_c + dy_c ]
	    +             wx0*wy0*wz1*values_c[i_c               + dz_c ]
	    +             wx1*wy0*wz1*values_c[i_c + dx_c        + dz_c ]
	    +             wx0*wy1*wz1*values_c[i_c        + dy_c + dz_c ]
	    +             wx1*wy1*wz1*values_c[i_c + dx_c + dy_c + dz_c ];
	}
      }
    }

    return (sizeof(T) * n3_c[0]*n3_c[1]*n3_c[2]);

  }
}

//======================================================================

