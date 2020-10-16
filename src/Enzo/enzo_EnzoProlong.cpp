// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlong.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of Enzo's prolongation

#include "enzo.hpp"

#include "cello.hpp"
// #define TRACE_PROLONG
// DEBUG_USE_LINEAR
#ifdef TRACE_PROLONG
#  undef TRACE_PROLONG
#  define TRACE_PROLONG(MSG) CkPrintf ("TRACE_PROLONG %s:%d %s\n",__FILE__,__LINE__,MSG);
#else
#  undef TRACE_PROLONG
#  define TRACE_PROLONG(MSG) /* ... */
#endif

//----------------------------------------------------------------------

EnzoProlong::EnzoProlong(std::string type,int positive) throw()
  : Prolong (),
    method_(-1),
    positive_(positive)
{
  TRACE_PROLONG("EnzoProlong()");
  if      (type == "3A")  method_ = 0;
  else if (type == "2A") method_ = 1;
  else if (type == "2B") method_ = 2;
  else {
    ERROR1("EnzoProlong::EnzoProlong",
	  "Unrecognized interpolation method %s",
	   type.c_str());
  }
}
//----------------------------------------------------------------------

void EnzoProlong::pup (PUP::er &p)
{
  TRACEPUP;

  Prolong::pup(p);

  p | method_;
}

//----------------------------------------------------------------------

int EnzoProlong::apply 
( precision_type precision,
  void *       values_f, int m3_f[3], int o3_f[3], int n3_f[3],
  const void * values_c, int m3_c[3], int o3_c[3], int n3_c[3],
  bool accumulate)
{
  TRACE_PROLONG("EnzoProlong::apply()");
#ifdef DEBUG_USE_LINEAR
  ProlongLinear p;
  return (p.apply(precision,values_f,m3_f,o3_f,n3_f,values_c,m3_c,o3_c,n3_c));
#else  
  return apply_((enzo_float *)     values_f,m3_f,o3_f,n3_f,
                (const enzo_float*)values_c,m3_c,o3_c,n3_c,accumulate);
#endif
}

//----------------------------------------------------------------------

int EnzoProlong::apply_
( 
 enzo_float * values_f, int m3_f[3], int o3_f[3], int n3_f[3],
 const enzo_float * values_c, int m3_c[3], int o3_c[3], int n3_c[3],
 bool accumulate)
{

  TRACE_PROLONG("EnzoProlong::apply_()");

  int rank = cello::rank();

  int error = 0;

  const int nxc=n3_c[0];
  const int nyc=n3_c[1];
  const int nzc=n3_c[2];
  const int mxc=m3_c[0];
  const int myc=m3_c[1];
  const int mzc=m3_c[2];
  const int mxbc=(rank>=1) ? nxc+2 : 1;
  const int mybc=(rank>=2) ? nyc+2 : 1;
  const int mzbc=(rank>=3) ? nzc+2 : 1;
  enzo_float *values_bc  = new enzo_float [mxbc*mybc*mzbc];
  if (rank == 2) {
    // copy array
    int ixc,iyc;
    int ixbc,iybc;
    const int ic0 = o3_c[0] + m3_c[0]*o3_c[1];
    for (iybc=0; iybc<mybc; iybc++) {
      iyc = iybc-1;
      iyc = std::max(iyc,0);
      iyc = std::min(iyc,nyc-1);
      for (ixbc=0; ixbc<mxbc; ixbc++) {
        ixc = ixbc-1;
        ixc = std::max(ixc,0);
        ixc = std::min(ixc,nxc-1);
        const int ibc=ixbc + mxbc*iybc;
        const int ic =ixc  +  mxc*iyc;
        values_bc[ibc] = values_c[ic+ic0];
      }
    }
  } else if (rank == 3) {
    ERROR("EnzoProlong","3D not implemented yet");
  }
  int ndim;
  enzo_float * parent;
  int pdims[3],pstart[3],pend[3];
  int r3[3] = {2,2,2};
  int gdims[3],gstart[3];

  parent = values_bc;

  pdims[0] = mxbc;
  pdims[1] = mybc;
  pdims[2] = mzbc;
  pstart[0] = 2;
  pstart[1] = 2;
  pstart[2] = 2;
  pend[0] = 2*nxc+1;
  pend[1] = 2*nyc+1;
  pend[2] = 2*nzc+1;

  const int nf = n3_f[0]*n3_f[1]*n3_f[2];
  enzo_float * grid = new enzo_float[nf];
  for (int i=0; i<nf; i++) grid[i] = -99.0;

  gdims[0] = n3_f[0];
  gdims[1] = n3_f[1];
  gdims[2] = n3_f[2];
  gstart[0] = 0;
  gstart[1] = 0;
  gstart[2] = 0;

  positive_=0;

  int size=gdims[0]/2 + 1;
  if (rank >= 2) size*=gdims[1]/2 + 1;
  if (rank >= 3) size*=gdims[2]/2 + 1;
  enzo_float * work = new enzo_float[size];

  FORTRAN_NAME(interpolate)
    (&rank,
     (enzo_float*)parent, pdims, pstart, pend, r3,
     grid, gdims, gstart, work, &method_,
     &positive_, &error);


  const int i3_f = o3_f[0] + m3_f[0]*(o3_f[1] + m3_f[1]*o3_f[2]);
  for (int iz=0; iz<n3_f[2]; iz++) {
    for (int iy=0; iy<n3_f[1]; iy++) {
      for (int ix=0; ix<n3_f[0]; ix++) {
        int iv = i3_f + ix+ m3_f[0]*(iy+ m3_f[1]*iz);
        int ig = ix+ n3_f[0]*(iy+ n3_f[1]*iz);
        values_f[iv] = grid[ig];
      }
    }
  }

  delete [] grid;
  delete [] work;
  delete [] values_bc;
  return (sizeof(enzo_float) * mxc*myc*mzc);
}  
