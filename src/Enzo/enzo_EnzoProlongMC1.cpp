// See LICENSE_CELLO file for license and copyright information

/// @file     field_EnzoProlongMC1.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of Enzo's prolongation

#include "enzo.hpp"

#include "cello.hpp"

EnzoProlongMC1::EnzoProlongMC1(std::string prolong_type) throw()
  : Prolong ()
{
  if      (prolong_type == "ThirdOrderA")  method_ = 0;
  else if (prolong_type == "SecondOrderA") method_ = 1;
  else if (prolong_type == "SecondOrderB") method_ = 2;
  else if (prolong_type == "SecondOrderC") method_ = 3;
  else if (prolong_type == "FirstOrderA")  method_ = 4;
  else {
    ERROR1("EnzoProlongMC1::EnzoProlongMC1",
	  "Unrecognized interpolation method %s",
	   prolong_type.c_str());
  }
}
//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void EnzoProlongMC1::pup (PUP::er &p)
{
  TRACEPUP;

  Prolong::pup(p);

  p | method_;
}

#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

int EnzoProlongMC1::apply 
( precision_type precision,
  void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
  const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3])
{
  switch (precision)  {

  case precision_single:

    return apply_( (float *)       values_f, nd3_f, im3_f, n3_f,
		   (const float *) values_c, nd3_c, im3_c, n3_c);

    break;

  case precision_double:

    return apply_( (double *)       values_f, nd3_f, im3_f, n3_f,
		   (const double *) values_c, nd3_c, im3_c, n3_c);

    break;

  default:

    ERROR1 ("EnzoProlongMC1::apply()",
	    "Unknown precision %d",
	    precision);

    return 0;
  }
}

//----------------------------------------------------------------------

template <class T>
int EnzoProlongMC1::apply_
( T *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
  const T * values_c, int nd3_c[3], int im3_c[3], int n3_c[3])
{
  int dx_c = 1;
  int dy_c = nd3_c[0];
  int dz_c = nd3_c[1];

  const double c1[4] = { 5.0*0.25, 3.0*0.25, 1.0*0.25, -1.0*0.25};
  const double c2[4] = {-1.0*0.25, 1.0*0.25, 3.0*0.25,  5.0*0.25};


  int rank = (nd3_f[2] > 1) ? 3 : ( (nd3_f[1] > 1) ? 2 : 1 );

  for (int i=0; i<rank; i++) {
    const char * xyz = "xyz";
    ASSERT3 ("EnzoProlongMC1::apply_",
	     "fine array %c-axis %d must be twice the size of the coarse axis %d",
	     xyz[i],n3_c[i],n3_f[i],
	     n3_f[i]==n3_c[i]*2);

    ASSERT2 ("EnzoProlongMC1::apply_",
	     "fine grid %c-axis %d must be divisible by 4",
	     xyz[i],n3_f[i],n3_f[i] % 4 == 0);
  }

  if (n3_f[1]==1) {

    for (int ix0=0; ix0<n3_f[0]; ix0+=4) {
      int ic_x = ix0/2;
      for (int ix=ix0; ix<ix0+4; ix++) {
	int icx = ix-ix0;
	int if_x = ix;

	int i_c = im3_c[0]+ic_x;
	int i_f = im3_f[0]+if_x;

	values_f[i_f] = ( c1[icx]*values_c[i_c] + 
			  c2[icx]*values_c[i_c+dx_c]);
      }
    }

    return (sizeof(T) * n3_c[0]);


  } else if (n3_f[2] == 1) {

    T min=100, avg = 0, sum = 0, max=-100;
    int ixmin=100,ixmax=-100,iymin=100,iymax=-100;
    int count = 0;
    for (int ix0=0; ix0<n3_f[0]; ix0+=4) {
      int ic_x = ix0/2;
      for (int iy0=0; iy0<n3_f[1]; iy0+=4) {
	int ic_y = iy0/2;

	for (int ix=ix0; ix<ix0+4; ix++) {
	  int icx = ix-ix0;
	  int if_x = ix;
	  for (int iy=iy0; iy<iy0+4; iy++) {
	    int icy = iy-iy0;
	    int if_y = iy;

	    int i_c = (im3_c[0]+ic_x) + nd3_c[0]*
	      (       (im3_c[1]+ic_y));
	    int i_f = (im3_f[0]+if_x) + nd3_f[0]*
	      (       (im3_f[1]+if_y));

	    ixmax=std::max(ixmax,im3_f[0]+if_x);
	    iymax=std::max(iymax,im3_f[1]+if_y);
	    ixmin=std::min(ixmin,im3_f[0]+if_x);
	    iymin=std::min(iymin,im3_f[1]+if_y);

	    values_f[i_f] = 
	      ( c1[icx]*c1[icy]*values_c[i_c] +
		c2[icx]*c1[icy]*values_c[i_c+dx_c] +
		c1[icx]*c2[icy]*values_c[i_c     +dy_c] +
		c2[icx]*c2[icy]*values_c[i_c+dx_c+dy_c]);
	    min=std::min(min,values_f[i_f]);
	    sum += values_f[i_f];
	    max=std::max(max,values_f[i_f]);
	    count++;
	    if (values_f[i_f] < 0.0) {
	      printf ("EnzoProlongMC1 i_c = %d %d\n",
		      (im3_c[0]+ic_x),(im3_c[1]+ic_y));
	      printf ("EnzoProlongMC1 i_f = %d %d\n",
		      (im3_f[0]+if_x),(im3_f[1]+if_y));
	      printf ("EnzoProlongMC1 icx,icy = %d %d\n",icx,icy);
	      printf ("EnzoProlongMC1 c1 = %f %f\n",c1[icx],c1[icy]);
	      printf ("EnzoProlongMC1 c2 = %f %f\n",c2[icx],c2[icy]);
	      printf ("EnzoProlongMC1 values_c = %f %f %f %f\n",
		      values_c[i_c],
		      values_c[i_c+dx_c],
		      values_c[i_c+dy_c],
		      values_c[i_c+dx_c+dy_c]);
	      printf ("EnzoProlongMC1 values_f = %f\n",values_f[i_f]);
	      
	    }
	  }
	}
      }
    }
    return (sizeof(T) * n3_c[0]*n3_c[1]);

  } else {

    for (int ix0=0; ix0<n3_f[0]; ix0+=4) {
      int ic_x = ix0/2;
      for (int iy0=0; iy0<n3_f[1]; iy0+=4) {
	int ic_y = iy0/2;
	for (int iz0=0; iz0<n3_f[2]; iz0+=4) {
	  int ic_z = iz0/2;

	  for (int ix=ix0; ix<ix0+4; ix++) {
	    int icx = ix-ix0;
	    int if_x = ix;
	    for (int iy=iy0; iy<iy0+4; iy++) {
	      int icy = iy-iy0;
	      int if_y = iy;
	      for (int iz=iz0; iz<iz0+4; iz++) {
		int icz = iz-iz0;
		int if_z = iz;

		int i_c = (im3_c[0]+ic_x) + nd3_c[0]*
		  (       (im3_c[1]+ic_y) + nd3_c[1]*
			  (im3_c[2]+ic_z));
		int i_f = (im3_f[0]+if_x) + nd3_f[0]*
		  (       (im3_f[1]+if_y) + nd3_f[1]*
			  (im3_f[2]+if_z));

		values_f[i_f] = 
		  ( c1[icx]*c1[icy]*c1[icz]*values_c[i_c] +
		    c2[icx]*c1[icy]*c1[icz]*values_c[i_c+dx_c] +
		    c1[icx]*c2[icy]*c1[icz]*values_c[i_c     +dy_c] +
		    c2[icx]*c2[icy]*c1[icz]*values_c[i_c+dx_c+dy_c] +
		    c1[icx]*c1[icy]*c2[icz]*values_c[i_c          +dz_c] +
		    c2[icx]*c1[icy]*c2[icz]*values_c[i_c+dx_c     +dz_c] +
		    c1[icx]*c2[icy]*c2[icz]*values_c[i_c     +dy_c+dz_c] +
		    c2[icx]*c2[icy]*c2[icz]*values_c[i_c+dx_c+dy_c+dz_c]);
	      }
	    }
	  }
	}
      }
    }
  }
  return (sizeof(T) * n3_c[0]*n3_c[1]*n3_c[2]);

}  
