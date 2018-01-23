// See LICENSE_CELLO file for license and copyright information

/// @file     test_ProlongLinear.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2018-01-19
/// @brief    Test program for the ProlongLinear class

#include "main.hpp"
#include "test.hpp"
#include <math.h>
#include "mesh.hpp"

double fun (double x, double y=0.0, double z=0.0)
{
  return 13.25
    - 7.5*x + 5.75*y - 1.125*z
    + 3.5*x*y + 1.875*y*z - 2.0*z*x
    + 7.125*x*y*z;
}

//  ghost = 0                ghost = 1
//
// 0   1   2   3            0   1   2   3
// o   o   o   o            o   o   o   o
//
//   O       O         O      O       O   
//   0       1
//  0.5     2.5      -1.5    0.5     2.5
//
//  p_c(i) - 2*g


double p_c(int i) { return 2.0*i + 0.5; }
double p_f(int i) { return i; }

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ProlongLinear");

  ProlongLinear * prolong = new ProlongLinear;

  unit_assert (prolong != NULL);

  int m3_f[3],n3_f[3],i3_f[3];
  int m3_c[3],n3_c[3],i3_c[3];

  double * v_c;
  double * v_f;
  //--------------------------------------------------


  for (int gx=0; gx<2; gx++) {

    char buffer[40+1];
    snprintf (buffer,40,"apply() 1D (%d)",gx);
    unit_func (buffer);

    m3_f[0] = 26;  m3_c[0] = 16;
    m3_f[1] = 1;   m3_c[1] = 1;
    m3_f[2] = 1;   m3_c[2] = 1;

    i3_f[0] = 3;  i3_c[0] = 2-gx;
    i3_f[1] = 0;  i3_c[1] = 0;
    i3_f[2] = 0;  i3_c[2] = 0;

    n3_f[0] = 20; n3_c[0] = 10+2*gx;
    n3_f[1] = 1;  n3_c[1] = 1;
    n3_f[2] = 1;  n3_c[2] = 1;

    v_c = new double [m3_c[0]*m3_c[1]*m3_c[2]];
    v_f = new double [m3_f[0]*m3_f[1]*m3_f[2]];

    std::fill_n(v_c,m3_c[0]*m3_c[1]*m3_c[2],0.0);
    std::fill_n(v_f,m3_f[0]*m3_f[1]*m3_f[2],0.0);
  
    for (int ix_c=i3_c[0]; ix_c<n3_c[0]+i3_c[0]; ix_c++) {
      int i=ix_c;
      double x = p_c(ix_c-i3_c[0]-gx);
      v_c[i] = fun(x);
    }

    prolong->apply (precision_double,
		    v_f, m3_f, i3_f, n3_f,
		    v_c, m3_c, i3_c, n3_c);

    for (int ix_f=i3_f[0]; ix_f<n3_f[0]+i3_f[0]; ix_f++) {
      int i=ix_f;
      double x = p_f(ix_f-i3_f[0]);
      unit_assert (v_f[i] == fun(x));
    }

    delete v_c;
    delete v_f;

  }

  //--------------------------------------------------

  unit_func ("apply() 2D ghost == 0");

  for (int gy=0; gy<2; gy++) {
    for (int gx=0; gx<2; gx++) {
      char buffer[40+1];
      snprintf (buffer,40,"apply() 2D (%d,%d)",gx,gy);
      unit_func (buffer);
    
      m3_f[0] = 26;  m3_c[0] = 16;
      m3_f[1] = 32;  m3_c[1] = 18;
      m3_f[2] = 1;   m3_c[2] = 1;

      i3_f[0] = 3;  i3_c[0] = 3-gx;
      i3_f[1] = 2;  i3_c[1] = 1-gy;
      i3_f[2] = 0;  i3_c[2] = 0;

      n3_f[0] = 20; n3_c[0] = 10+2*gx;
      n3_f[1] = 28; n3_c[1] = 14+2*gy;
      n3_f[2] = 1;  n3_c[2] = 1;

      v_c = new double [m3_c[0]*m3_c[1]*m3_c[2]];
      v_f = new double [m3_f[0]*m3_f[1]*m3_f[2]];

      std::fill_n(v_c,m3_c[0]*m3_c[1]*m3_c[2],0.0);
      std::fill_n(v_f,m3_f[0]*m3_f[1]*m3_f[2],0.0);
  
      for (int iy_c=i3_c[1]; iy_c<n3_c[1]+i3_c[1]; iy_c++) {
	double y = p_c(iy_c-i3_c[1]-gy);
	for (int ix_c=i3_c[0]; ix_c<n3_c[0]+i3_c[0]; ix_c++) {
	  int i=ix_c + m3_c[0]*iy_c;
	  double x = p_c(ix_c-i3_c[0]-gx);
	  v_c[i] = fun(x,y);
	}
      }

      prolong->apply (precision_double,
		      v_f, m3_f, i3_f, n3_f,
		      v_c, m3_c, i3_c, n3_c);
  
      for (int iy_f=i3_f[1]; iy_f<n3_f[1]+i3_f[1]; iy_f++) {
	double y = p_f(iy_f-i3_f[1]);
	for (int ix_f=i3_f[0]; ix_f<n3_f[0]+i3_f[0]; ix_f++) {
	  int i=ix_f + m3_f[0]*iy_f;
	  double x = p_f(ix_f-i3_f[0]);
	  unit_assert (v_f[i] == fun(x,y));
	}
      }

      delete v_c;
      delete v_f;
    }
  }

  //--------------------------------------------------

  for (int gz=0; gz<2; gz++) {
    for (int gy=0; gy<2; gy++) {
      for (int gx=0; gx<2; gx++) {
	char buffer[40+1];
	snprintf (buffer,40,"apply() 3D (%d,%d,%d)",gx,gy,gz);
	unit_func (buffer);
    
	m3_f[0] = 16;  m3_c[0] = 13;
	m3_f[1] = 19;  m3_c[1] = 10;
	m3_f[2] = 17;  m3_c[2] = 30;

	i3_f[0] = 3;  i3_c[0] = 3-gx;
	i3_f[1] = 2;  i3_c[1] = 1-gy;
	i3_f[2] = 1;  i3_c[2] = 2-gz;

	n3_f[0] = 10; n3_c[0] = 5+2*gx;
	n3_f[1] = 12; n3_c[1] = 6+2*gy;
	n3_f[2] = 12; n3_c[2] = 6+2*gz;

	v_c = new double [m3_c[0]*m3_c[1]*m3_c[2]];
	v_f = new double [m3_f[0]*m3_f[1]*m3_f[2]];

	std::fill_n(v_c,m3_c[0]*m3_c[1]*m3_c[2],0.0);
	std::fill_n(v_f,m3_f[0]*m3_f[1]*m3_f[2],0.0);
  
	for (int iz_c=i3_c[2]; iz_c<n3_c[2]+i3_c[2]; iz_c++) {
	  double z = p_c(iz_c-i3_c[2]-gz);
	  for (int iy_c=i3_c[1]; iy_c<n3_c[1]+i3_c[1]; iy_c++) {
	    double y = p_c(iy_c-i3_c[1]-gy);
	    for (int ix_c=i3_c[0]; ix_c<n3_c[0]+i3_c[0]; ix_c++) {
	      int i=ix_c + m3_c[0]*(iy_c + m3_c[1]*iz_c);
	      double x = p_c(ix_c-i3_c[0]-gx);
	      v_c[i] = fun(x,y,z);
	    }
	  }
	}

	prolong->apply (precision_double,
			v_f, m3_f, i3_f, n3_f,
			v_c, m3_c, i3_c, n3_c);
  
	for (int iz_f=i3_f[2]; iz_f<n3_f[2]+i3_f[2]; iz_f++) {
	  double z = p_f(iz_f-i3_f[2]);
	  for (int iy_f=i3_f[1]; iy_f<n3_f[1]+i3_f[1]; iy_f++) {
	    double y = p_f(iy_f-i3_f[1]);
	    for (int ix_f=i3_f[0]; ix_f<n3_f[0]+i3_f[0]; ix_f++) {
	      int i=ix_f + m3_f[0]*(iy_f + m3_f[1]*iz_f);
	      double x = p_f(ix_f-i3_f[0]);
	      unit_assert (v_f[i] == fun(x,y,z));
	    }
	  }
	}

	delete v_c;
	delete v_f;
      }
    }
  }

  //--------------------------------------------------
  
  delete prolong;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

