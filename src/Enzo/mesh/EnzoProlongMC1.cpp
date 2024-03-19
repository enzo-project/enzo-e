// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlongMC1.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of Enzo's prolongation (INCOMPLETE)

#include "Enzo/mesh/mesh.hpp"
#include "Enzo/enzo.hpp"
#include "Cello/cello.hpp"

EnzoProlongMC1::EnzoProlongMC1(std::string type) throw()
  : Prolong (),
    method_(-1)
{
  if      (type == "ThirdOrderA")  method_ = 0;
  else if (type == "SecondOrderA") method_ = 1;
  else if (type == "SecondOrderB") method_ = 2;
  else if (type == "SecondOrderC") method_ = 3;
  else if (type == "FirstOrderA")  method_ = 4;
  else {
    ERROR1("EnzoProlongMC1::EnzoProlongMC1",
	  "Unrecognized interpolation method %s",
	   type.c_str());
  }
}
//----------------------------------------------------------------------

void EnzoProlongMC1::pup (PUP::er &p)
{
  TRACEPUP;

  Prolong::pup(p);

  p | method_;
}

//----------------------------------------------------------------------

void EnzoProlongMC1::apply 
( precision_type precision,
  void *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
  const void * values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
  bool accumulate)
{
  apply_( (enzo_float * )       values_f, nd3_f, im3_f, n3_f,
          (const enzo_float * ) values_c, nd3_c, im3_c, n3_c,
          accumulate);
}

//----------------------------------------------------------------------

void EnzoProlongMC1::apply_
( enzo_float *       values_f, int nd3_f[3], int im3_f[3], int n3_f[3],
  const enzo_float * values_c, int nd3_c[3], int im3_c[3], int n3_c[3],
  bool accumulate)
{
  ASSERT ("EnzoProlongPoisson::apply_",
	  "accumulate=true is not implemented yet",
	  ! accumulate);
  
  const int dx_c = 1;
  const int dy_c = nd3_c[0];

  const int dx_f = 1;
  const int dy_f = nd3_f[0];

  const int rank = (nd3_f[2] > 1) ? 3 : ( (nd3_f[1] > 1) ? 2 : 1 );

  const int ndc = nd3_c[0]*nd3_c[1]*nd3_c[2];

  enzo_float * work = new enzo_float [ ndc ];

  for (int i=0; i<rank; i++) {
    const char * xyz = "xyz";
    ASSERT3 ("EnzoProlongMC1::apply_",
	     "fine array %c-axis %d must be twice the size of the coarse axis %d",
	     xyz[i],n3_c[i],n3_f[i],
	     n3_f[i]==n3_c[i]*2);
  }

  if (n3_f[1]==1) {

    for (int ix=0; ix<n3_c[0]-1; ix++) {
      work[ix] = 0.5*(values_c[ix] + values_c[ix+dx_c]);
    }

    delete [] work;

  } else if (n3_f[2] == 1) {

    for (int ix_c=0; ix_c<n3_c[0]; ix_c++) {
      for (int iy_c=0; iy_c<n3_c[1]; iy_c++) {
	int i_c = ix_c + n3_c[0]*iy_c;
	work[i_c] = 0.25*(values_c[i_c] + 
			  values_c[i_c+dx_c] +
			  values_c[i_c+dy_c] +
			  values_c[i_c+dx_c+dy_c]);
      }
    }

    for (int ix_c=0; ix_c<n3_c[0]; ix_c++) {
      int ix_f = 2*ix_c;
      for (int iy_c=0; iy_c<n3_c[1]; iy_c++) {
	int iy_f = 2*iy_c;
	int i_c = ix_c + n3_c[0]*iy_c;
	int i_f = ix_f + n3_f[0]*iy_f;

	int i00_c = i_c;
	int i10_c = i_c + dx_c;
	int i01_c = i_c +        dy_c;
	int i11_c = i_c + dx_c + dy_c;

	enzo_float fb = values_c[i00_c];

	enzo_float df0 = std::min(fabs(fb - work[i00_c]),fabs(fb - work[i11_c]));
	enzo_float df1 = std::min(fabs(fb - work[i10_c]),fabs(fb - work[i01_c]));

	df0 = (fb-work[i00_c]) < 0.0 ? -df0 : df0;
	df1 = (fb-work[i10_c]) < 0.0 ? -df1 : df1;

	if ((work[i11_c]-fb)*(fb-work[i00_c]) < 0.0) df0 = 0.0;
	if ((work[i01_c]-fb)*(fb-work[i10_c]) < 0.0) df1 = 0.0;

	int i00_f = i_f;
	int i10_f = i_f + dx_f;
	int i01_f = i_f +        dy_f;
	int i11_f = i_f + dx_f + dy_f;

	enzo_float fx = df0 - df1;
	enzo_float fy = df0 + df1;

	values_f[i00_f]	= fb - 0.5*fx - 0.5*fy;
	values_f[i10_f]	= fb + 0.5*fx - 0.5*fy;
	values_f[i01_f]	= fb - 0.5*fx + 0.5*fy;
	values_f[i11_f]	= fb + 0.5*fx + 0.5*fy;

      }
    }

    delete [] work;
  
  } else {

    ERROR("EnzoProlongMC1",
	  "3D not implemented yet");

    delete [] work;
  }

}  
