// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlong.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of Enzo's prolongation

#include "enzo.hpp"

#include "cello.hpp"

// #define TRACE_PROLONG
// #define DEBUG_ENZO_PROLONG

// #define DEBUG_ARRAY

#ifdef TRACE_PROLONG
#  undef TRACE_PROLONG
#  define TRACE_PROLONG(MSG) CkPrintf ("TRACE_PROLONG %s:%d %s\n",__FILE__,__LINE__,MSG);
#else
#  undef TRACE_PROLONG
#  define TRACE_PROLONG(MSG) /* ... */
#endif

#ifdef DEBUG_ARRAY

#   define DEBUG_PRINT_ARRAY0(NAME,ARRAY,m3,n3,o3)              \
  DEBUG_PRINT_ARRAY_(NAME,ARRAY,m3,n3,o3[0],o3[1],o3[2])
#   define DEBUG_PRINT_ARRAY(NAME,ARRAY,m3,n3)  \
  DEBUG_PRINT_ARRAY_(NAME,ARRAY,m3,n3,0,0,0)
#   define DEBUG_PRINT_ARRAY_(NAME,ARRAY,m3,n3,ox,oy,oz)        \
  {                                                             \
    CkPrintf ("PADDED_ARRAY_VALUES %s:%d %s %p\n",              \
              __FILE__,__LINE__,NAME,(void*)ARRAY);             \
    CkPrintf ("m3 %d %d %d n3 %d %d %d o3 %d %d %d\n",          \
              m3[0],m3[1],m3[2],n3[0],n3[1],n3[2],ox,oy,oz); \
    const int o = ox + m3[0]*(oy + m3[1]*oz);                   \
    for (int iz=0; iz<n3[2]; iz++) {                            \
      for (int iy=0; iy<n3[1]; iy++) {                          \
        CkPrintf ("PADDED_ARRAY_VALUES %s %p %d %d %d: ",       \
                  NAME,(void*)ARRAY,0,iy,iz);                   \
        for (int ix=0; ix<n3[0]; ix++) {                        \
          int i = ix+ m3[0]*(iy+ m3[1]*iz);                     \
          CkPrintf (" %6.3g",ARRAY[o+i]);                       \
        }                                                       \
        CkPrintf ("\n");                                        \
      }                                                         \
    }                                                           \
  }

#   define DEBUG_FILL_ARRAY0(NAME,ARRAY,m3,n3,o3)       \
  DEBUG_FILL_ARRAY_(NAME,ARRAY,m3,n3,o3[0],o3[1],o3[2])
#   define DEBUG_FILL_ARRAY(NAME,ARRAY,m3,n3)   \
  DEBUG_FILL_ARRAY_(NAME,ARRAY,m3,n3,0,0,0)
#   define DEBUG_FILL_ARRAY_(NAME,ARRAY,m3,n3,ox,oy,oz) \
  {                                                     \
    const int o = ox + m3[0]*(oy + m3[1]*oz);           \
    for (int iz=0; iz<n3[2]; iz++) {                    \
      for (int iy=0; iy<n3[1]; iy++) {                  \
        for (int ix=0; ix<n3[0]; ix++) {                \
          int i = ix+ m3[0]*(iy+ m3[1]*iz);             \
          ARRAY[o+i] = o+i;                             \
        }                                               \
      }                                                 \
    }                                                   \
  }

#   define DEBUG_COPY_ARRAY0(NAME,ARRAY_D,ARRAY_S,m3,n3,o3)             \
  DEBUG_COPY_ARRAY_(NAME,ARRAY_D,ARRAY_S,m3,n3,o3[0],o3[1],o3[2])
#   define DEBUG_COPY_ARRAY(NAME,ARRAY_D,ARRAY_S,m3,n3) \
  DEBUG_COPY_ARRAY_(NAME,ARRAY_D,ARRAY_S,m3,n3,0,0,0)
#   define DEBUG_COPY_ARRAY_(NAME,ARRAY_D,ARRAY_S,m3,n3,ox,oy,oz)       \
  {                                                                     \
    const int o = ox + m3[0]*(oy + m3[1]*oz);                           \
    for (int iz=0; iz<n3[2]; iz++) {                                    \
      for (int iy=0; iy<n3[1]; iy++) {                                  \
        for (int ix=0; ix<n3[0]; ix++) {                                \
          int i = ix+ m3[0]*(iy+ m3[1]*iz);                             \
          ARRAY_D[o+i] = ARRAY_S[o+i];                                  \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }

#else

#   define DEBUG_PRINT_ARRAY0(NAME,ARRAY,m3,n3,o3)  /* ... */
#   define DEBUG_PRINT_ARRAY(NAME,ARRAY,m3,n3)  /* ... */
#   define DEBUG_PRINT_ARRAY_(NAME,ARRAY,m3,n3,ox,oy,oz)  /* ... */

#   define DEBUG_FILL_ARRAY0(NAME,ARRAY,m3,n3,o3)  /* ... */
#   define DEBUG_FILL_ARRAY(NAME,ARRAY,m3,n3)  /* ... */
#   define DEBUG_FILL_ARRAY_(NAME,ARRAY,m3,n3,ox,oy,oz)  /* ... */

#   define DEBUG_COPY_ARRAY0(NAME,ARRAY_D,ARRAY_S,m3,n3,o3)  /* ... */
#   define DEBUG_COPY_ARRAY(NAME,ARRAY_D,ARRAY_S,m3,n3)  /* ... */
#   define DEBUG_COPY_ARRAY_(NAME,ARRAY_D,ARRAY_S,m3,n3,ox,oy,oz)  /* ... */
#endif
//----------------------------------------------------------------------

EnzoProlong::EnzoProlong(std::string type,int positive) throw()
  : Prolong (),
    method_(-1),
    positive_(positive)
{
  TRACE_PROLONG("EnzoProlong()");
  if      (type == "3A")  method_ = 0;
  else if (type == "linear") method_ = -1;
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

void EnzoProlong::apply 
( precision_type precision,
  void *       values_f, int m3_f[3], int o3_f[3], int n3_f[3],
  const void * values_c, int m3_c[3], int o3_c[3], int n3_c[3],
  bool accumulate)
{
  if (method_ >= 0) {
    apply_((enzo_float *)     values_f,m3_f,o3_f,n3_f,
           (const enzo_float*)values_c,m3_c,o3_c,n3_c,accumulate);
  } else {
    WARNING("EnzoProlong::apply()","CALLING ProlongLinear");
    ProlongLinear prolong_linear;
    prolong_linear.apply
      (precision,
       values_f,m3_f,o3_f,n3_f,
       values_c,m3_c,o3_c,n3_c,accumulate);
  }
}
//----------------------------------------------------------------------

void EnzoProlong::apply_
( 
 enzo_float * values_f, int m3_f[3], int o3_f[3], int n3_f[3],
 const enzo_float * values_c, int m3_c[3], int o3_c[3], int n3_c[3],
 bool accumulate)
{
  int rank = cello::rank();

  int error = 0;

  int r3[3] = {2,2,2};
  int gdims[3],gstart[3],pdims[3],pstart[3],pend[3];

  //  enzo_float * grid = new enzo_float[nf];
  //  DEBUG_FILL_ARRAY("grid",grid,n3_f,n3_f);

  positive_ = 0;
  for (int i=0; i<rank; i++) {
    pdims[i] = m3_c[i];
    pstart[i] = 2;
    pend[i] = n3_f[i]+1;
    gdims[i] = m3_f[i];
    gstart[i] = 0;
  }
  
  int size=gdims[0]/2 + 1;
  if (rank >= 2) size*=gdims[1]/2 + 1;
  if (rank >= 3) size*=gdims[2]/2 + 1;
  enzo_float * work = new enzo_float[2*size+100];

#ifdef DEBUG_ENZO_PROLONG
  CkPrintf ("DEBUG_ENZO_PROLONG EnzoProlong\n");
  CkPrintf ("DEBUG_ENZO_PROLONG m3_f    %d %d %d\n",m3_f[0],m3_f[1],m3_f[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG n3_f    %d %d %d\n",n3_f[0],n3_f[1],n3_f[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG o3_f    %d %d %d\n",o3_f[0],o3_f[1],o3_f[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG m3_c    %d %d %d\n",m3_c[0],m3_c[1],m3_c[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG n3_c    %d %d %d\n",n3_c[0],n3_c[1],n3_c[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG o3_c    %d %d %d\n",o3_c[0],o3_c[1],o3_c[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG\n");

  // convert to fortran ordering for interpolate parameters
  // for (int i=0; i<rank/2; i++) {
  //   std::swap(pdims[i],pdims[rank-i-1]);
  //   std::swap(pstart[i],pstart[rank-i-1]);
  //   std::swap(pend[i],pend[rank-i-1]);
  //   std::swap(gdims[i],gdims[rank-i-1]);
  //   std::swap(gstart[i],gstart[rank-i-1]);
  // }
  
  CkPrintf ("DEBUG_ENZO_PROLONG pdims   %d %d %d\n",pdims[0],pdims[1],pdims[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG pstart  %d %d %d\n",pstart[0],pstart[1],pstart[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG pend    %d %d %d\n",pend[0],pend[1],pend[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG r3      %d %d %d\n",r3[0],r3[1],r3[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG gdims   %d %d %d\n",gdims[0],gdims[1],gdims[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG gstart  %d %d %d\n",gstart[0],gstart[1],gstart[2]);
  CkPrintf ("DEBUG_ENZO_PROLONG rank    %d\n",rank);

#endif

  int o_c = o3_c[0] + m3_c[0]*(o3_c[1] + m3_c[1]*o3_c[2]);
  int o_f = o3_f[0] + m3_f[0]*(o3_f[1] + m3_f[1]*o3_f[2]);

  FORTRAN_NAME(interpolate)
    (&rank,
     ((enzo_float*)(values_c))+o_c, pdims, pstart, pend, r3,
     ((enzo_float*)(values_f))+o_f, gdims, gstart, work, &method_,
     &positive_, &error);

  DEBUG_PRINT_ARRAY0("values_c",values_c,m3_c,n3_c,o3_c);
  DEBUG_PRINT_ARRAY0("values_f",values_f,m3_f,n3_f,o3_f);

  delete [] work;
}
