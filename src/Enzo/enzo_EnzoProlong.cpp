// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlong.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of Enzo's prolongation

#include "enzo.hpp"

#include "cello.hpp"
#define TRACE_PROLONG
// #define DEBUG_ENZO_PROLONG
// #define TRACE_PADDED_ARRAY_VALUES

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
  return apply_((enzo_float *)     values_f,m3_f,o3_f,n3_f,
                (const enzo_float*)values_c,m3_c,o3_c,n3_c,accumulate);
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

  int r3[3] = {2,2,2};
  int gdims[3],gstart[3];


  int pdims[3] =
    { (rank>=1) ? n3_c[0] : 1,
      (rank>=2) ? n3_c[1] : 1,
      (rank>=3) ? n3_c[2] : 1};
  int pstart[3] = {2,2,2};
  
  int pend[3] = {(rank >= 1) ? 2*n3_c[0]-3 : 3,
                 (rank >= 2) ? 2*n3_c[1]-3 : 3,
                 (rank >= 3) ? 2*n3_c[2]-3 : 3};

  const int nf = n3_f[0]*n3_f[1]*n3_f[2];

  const enzo_float * parent(values_c);
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

#ifdef DEBUG_ENZO_PROLONG
  CkPrintf ("DEBUG_ENZO_PROLONG EnzoProlong\n");
  //  CkPrintf ("DEBUG_ENZO_PROLONG NEW n3_f    %d %d %d\n",n3_f[0],n3_f[1],n3_f[2]);
  //  CkPrintf ("DEBUG_ENZO_PROLONG NEW n3_c    %d %d %d\n",n3_c[0],n3_c[1],n3_c[2]);
  //  CkPrintf ("DEBUG_ENZO_PROLONG NEW o3_c    %d %d %d\n",o3_c[0],o3_c[1],o3_c[2]);

  CkPrintf ("DEBUG_ENZO_PROLONG NEW pdims   %d %d %d\n",pdims[0],pdims[1],pdims[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG NEW pstart  %d %d %d\n",pstart[0],pstart[1],pstart[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG NEW pend    %d %d %d\n",pend[0],pend[1],pend[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG NEW r3      %d %d %d\n",r3[0],r3[1],r3[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG NEW gdims   %d %d %d\n",gdims[0],gdims[1],gdims[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG NEW gstart  %d %d %d\n",gstart[0],gstart[1],gstart[2]);

  CkPrintf ("TRACE_PADDED_ARRAY EnzoProlong::apply parent %p sum %g\n",
            (void*)parent,
            cello::sum(parent,
                       m3_c[0],m3_c[1],m3_c[2],
                       o3_c[0],o3_c[1],o3_c[2],
                       n3_c[0],n3_c[1],n3_c[2]));

#endif

#ifdef DEBUG_ENZO_PROLONG
  {
    int im=o3_c[0]+m3_c[0]*(o3_c[1]+m3_c[1]*o3_c[2]);
    int ip=(o3_c[0]+n3_c[0]-1)+m3_c[0]*((o3_c[1]+n3_c[1]-1)+m3_c[1]*(o3_c[2]+n3_c[2]-1));
    
  CkPrintf ("DEBUG_PROLONG %s:%d parent first %g\n",
            __FILE__,__LINE__,parent[im]);
  CkPrintf ("DEBUG_PROLONG %s:%d parent sum %g\n",
            __FILE__,__LINE__,cello::sum
            (parent,
             m3_c[0],m3_c[1],m3_c[2],
             o3_c[0],o3_c[1],o3_c[2],
             n3_c[0],n3_c[1],n3_c[2]));
  CkPrintf ("DEBUG_PROLONG %s:%d parent last %g\n",
            __FILE__,__LINE__,parent[ip]);
  }

#endif  

#ifdef TRACE_PADDED_ARRAY_VALUES
  CkPrintf ("PADDED_ARRAY_VALUES %s:%d parent %p\n",
            __FILE__,__LINE__,(void*)parent);
  for (int iz=0; iz<n3_c[2]; iz++) {
    for (int iy=0; iy<n3_c[1]; iy++) {
      CkPrintf ("PADDED_ARRAY_VALUES parent %p %d %d %d: ",parent,0,iy,iz);
      for (int ix=0; ix<n3_c[0]; ix++) {
        int i = ix+ n3_c[0]*(iy+ n3_c[1]*iz);
        CkPrintf (" %6.3g",parent[i]);
      }
      CkPrintf ("\n");
    }
  }
#endif
  FORTRAN_NAME(interpolate)
    (&rank,
     (enzo_float*)parent, pdims, pstart, pend, r3,
     grid, gdims, gstart, work, &method_,
     &positive_, &error);


#ifdef TRACE_PADDED_ARRAY_VALUES
  CkPrintf ("PADDED_ARRAY_VALUES %s:%d grid %p\n",
            __FILE__,__LINE__,(void*)grid);
  for (int iz=0; iz<n3_f[2]; iz++) {
    for (int iy=0; iy<n3_f[1]; iy++) {
      CkPrintf ("PADDED_ARRAY_VALUES grid %p %d %d %d: ",grid,0,iy,iz);
      for (int ix=0; ix<n3_f[0]; ix++) {
        int i = ix+ n3_f[0]*(iy+ n3_f[1]*iz);
        CkPrintf (" %6.3g",grid[i]);
      }
      CkPrintf ("\n");
    }
  }
#endif

#ifdef DEBUG_ENZO_PROLONG
  {
    int im=o3_f[0]+m3_f[0]*(o3_f[1]+m3_f[1]*o3_f[2]);
    int ip=(o3_f[0]+n3_f[0]-1)+m3_f[0]*((o3_f[1]+n3_f[1]-1)+m3_f[1]*(o3_f[2]+n3_f[2]-1));
    CkPrintf ("DEBUG_PROLONG %s:%d grid first %p %d %g\n",
              __FILE__,__LINE__,(void*)grid,im,grid[im]);
    CkPrintf ("DEBUG_PROLONG %s:%d grid sum %g\n",
              __FILE__,__LINE__,cello::sum
              (grid,
               m3_c[0],m3_c[1],m3_c[2],
               o3_c[0],o3_c[1],o3_c[2],
               n3_c[0],n3_c[1],n3_c[2]));
    // VALGRIND INVALID READ 
    // CkPrintf ("DEBUG_PROLONG %s:%d grid last %g\n",
    //           __FILE__,__LINE__,grid[ip]);
  }
#endif  
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
  return (sizeof(enzo_float) * m3_c[0]*m3_c[1]*m3_c[2]);
  
}  
