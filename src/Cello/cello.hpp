// See LICENSE_CELLO file for license and copyright information

#ifndef CELLO_HPP
#define CELLO_HPP

/// @file    cello.hpp
/// @author  James Bordner (jobordner@ucsd.edu)
/// @date    Thu Nov 11 17:08:38 PST 2010
/// @brief   Include Cello global configuration settings
///
/// This file includes system includes; defines global template functions
/// such as MIN(), MAX(), and INDEX(); and global enumerated types.
/// It also initializes precision-related defines, including default_precision
/// Some global functions and constants are define with cello namespace,
/// sich as cello:pi and cello::err_rel(), etc.

//----------------------------------------------------------------------
// SYSTEM INCLUDES
//----------------------------------------------------------------------

// (check CMakeLists.txt to check to see if headers in the associated
//  precompiled header needs to change whenever any of these include statements
//  are removed - the precompiled header should only contain a subset of them)
#include <execinfo.h>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstdlib>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits> // std::remove_cv, std::is_same_v
#include <vector>

#include <charm++.h>

#include "pup_stl.h"

#include "cello_defines.hpp"
#include "cello_Sync.hpp"

// #define DEBUG_CHECK

class Block;
class Boundary;
class Config;
class CProxy_Block;
class Factory;
class Field;
class FieldDescr;
class Grouping;
class Hierarchy;
class Monitor;
class Output;
class Parameters;
class ParticleDescr;
class Problem;
class Refresh;
class ScalarDescr;
class Simulation;
class Solver;
class Units;

#ifdef CELLO_DEBUG
# define TRACE_ONCE                                             \
  {                                                             \
    static bool first = true;                                   \
    if (first) {                                                \
      first = false;                                            \
      CkPrintf ("TRACE_ONCE %s:%d\n",__FILE__,__LINE__);        \
    }                                                           \
  }
#else
# define TRACE_ONCE /* ... */
#endif

//----------------------------------------------------------------------
// TEMPLATE FUNCTIONS
//----------------------------------------------------------------------

template <class T>
inline T MIN(const T &a, const T &b) 
{  return a < b ? a : b; }

template <class T>
inline T MAX(const T &a, const T &b) 
{  return a > b ? a : b; }

inline int INDEX(int ix,int iy,int iz,int nx,int ny) 
{  return ix+nx*(iy+ny*iz); }

inline int INDEX2(int ix,int iy,int nx) 
{  return ix+nx*iy; }

template <class T>
inline void SWAP(T &a, T &b) 
{  T t = a; a=b; b=t; }

//----------------------------------------------------------------------
// ENUMERATED TYPES
//----------------------------------------------------------------------

/// @enum     face_enum
/// @brief    Face [lower|upper]
enum face_enum {
  face_lower = 0,
  face_upper = 1,
  face_all
};

/// @enum     axis_enum
/// @brief    Axis [x|y|z]
enum axis_enum {
  axis_x = 0,
  axis_y = 1,
  axis_z = 2,
  axis_all
};
typedef int axis_type;

/// @enum     refresh_type
/// @brief    Type of mesh refinement--coarsen, refine, or stay the same
enum refresh_type {
  refresh_unknown,
  refresh_coarse,
  refresh_same,
  refresh_fine
};

/// @enum     reduce_enum
/// @brief    Reduction operator, used for image projections
enum reduce_enum {
  reduce_unknown, /// Unknown reduction
  reduce_min,     /// Minimal value along the axis
  reduce_max,     /// Maximal value along the axis
  reduce_avg,     /// Average value along the axis
  reduce_sum,     /// Sum of values along the axis
  reduce_set      /// Value of last processed (used for mesh plotting)
};

typedef int reduce_type;


/// @enum precision_enum
/// @brief list of known floating-point precision, used for Field
enum precision_enum {
  // @@@ KEEP IN SYNCH WITH precision_name in cello.cpp
  precision_unknown,     //  unknown precision
  precision_default,     //  default precision
  precision_single,      //  32-bit field data
  precision_double,      //  64-bit field data
  precision_extended80,  //  80-bit field data
  precision_extended96,  //  96-bit field data
  precision_quadruple,   // 128-bit field data
};
typedef int precision_type;

#ifdef CONFIG_PRECISION_SINGLE
   typedef float       cello_float;
#elif  CONFIG_PRECISION_DOUBLE
   typedef double      cello_float;
#elif  CONFIG_PRECISION_QUAD
   typedef long double cello_float;
#else
#  error "Must define one of CONFIG_PRECISION_[SINGLE|DOUBLE|QUAD]"
#endif

/// @enum type_enum
/// @brief list of known scalar types, including ints as well as floats, used for Field types and Particle attributes
enum type_enum {
  type_unknown,     // unknown type
  type_default,    // "default" floating-point precision, e.g. enzo_float
  type_single,
  type_float = type_single,
  type_double,
  type_extended80,
  type_extended96,
  type_quadruple,
  type_long_double = type_quadruple,
  type_int8,
  type_char = type_int8,
  type_int16,
  type_short = type_int16,
  type_int32,
  type_int = type_int32,
  type_int64,
  type_long_long = type_int64,
  NUM_TYPES
};

#ifdef CONFIG_PRECISION_SINGLE
#   define default_precision precision_single
#   define default_type      type_single
#   define SCALAR_DEFINED
#   define default_precision_string "single"
#endif
#ifdef CONFIG_PRECISION_DOUBLE
#   ifdef SCALAR_DEFINED
#      define SCALAR_ERROR
#   endif
#   define default_precision precision_double
#   define default_type      type_double
#   define SCALAR_DEFINED
#   define default_precision_string "double"
#endif
#ifdef CONFIG_PRECISION_QUAD
#   ifdef SCALAR_DEFINED
#      define SCALAR_ERROR
#   endif
#   define default_precision precision_quad
#   define default_type      type_quad
#   define SCALAR_DEFINED
#   define default_precision_string "quadruple"
#endif

#ifndef SCALAR_DEFINED
#   error None of CONFIG_PRECISION_[SINGLE|DOUBLE|QUAD] defined
#endif

#ifdef ERROR_SCALAR
#   error Multiple CONFIG_PRECISION_[SINGLE|DOUBLE|QUAD] defined
#endif

enum class MsgType { msg_refine, msg_check };

/// Length of hex message tags used for debugging
#define TAG_LEN 8

/// represents a kind of initial cycle
///
/// @note
/// At the moment, we distinguish between a charm-based restart and a
/// non-charm-based restart for completeness. But to my knowledge, we can't
/// actually identify a charm-based restart
enum class InitCycleKind {
  fresh,                    ///< doesn't follow a restart
  charmrestart,             ///< follows a charm-based restart
  fresh_or_noncharm_restart ///< any kind other than a charm-based restart
};

//----------------------------------------------------------------------
/// Macros for debugging
//----------------------------------------------------------------------

#define TRACE_PROLONG_SUM(FUNC,V_F,M_F,O_F,N_F,V_C,M_C,O_C,N_C,A)       \
  {                                                                     \
    long double sum_f=0.0, sum_c=0.0;                                   \
    int count_f=0,count_c=0;                                            \
    int ic0 = O_C[0] + M_C[0]*(O_C[1] + M_C[1]*O_C[2]);                 \
    int if0 = O_F[0] + M_F[0]*(O_F[1] + M_F[1]*O_F[2]);                 \
    for (int icz=0; icz<N_C[2]; icz++) {                                \
      for (int icy=0; icy<N_C[1]; icy++) {                              \
        for (int icx=0; icx<N_C[0]; icx++) {                            \
          int i_c=ic0 + icx + M_C[0]*(icy + M_C[1]*icz);                \
          sum_c+=V_C[i_c];                                              \
          count_c++;                                                    \
        }                                                               \
      }                                                                 \
    }                                                                   \
    for (int ifz=0; ifz<N_F[2]; ifz++) {                                \
      for (int ify=0; ify<N_F[1]; ify++) {                              \
        for (int ifx=0; ifx<N_F[0]; ifx++) {                            \
          int i_f=if0 + ifx + M_F[0]*(ify + M_F[1]*ifz);                \
          sum_f+=V_F[i_f];                                              \
          count_f++;                                                    \
        }                                                               \
      }                                                                 \
    }                                                                   \
    CkPrintf ("TRACE_PROLONG %d %s\n",A,FUNC.c_str());                  \
    CkPrintf ("TRACE_PROLONG %d mc3 %d %d nc3 %d %d oc3 %d %d\n",       \
              A,M_C[0],M_C[1],N_C[0],N_C[1],O_C[0],O_C[1]);             \
    CkPrintf ("TRACE_PROLONG %d mf3 %d %d nf3 %d %d of3 %d %d\n",       \
              A,M_F[0],M_F[1],N_F[0],N_F[1],O_F[0],O_F[1]);             \
    CkPrintf ("TRACE_PROLONG %d sum_c %d %20.15Lg  sum_f %d %20.15Lg\n", \
              A,count_c,sum_c,count_f,sum_f);                           \
  }

//----------------------------------------------------------------------

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
/// Macros for sizing, saving, and restoring data from buffers
//----------------------------------------------------------------------

#define SIZE_SCALAR_TYPE(COUNT,TYPE,VALUE)      \
  {						\
    (COUNT) += sizeof(TYPE);			\
  }

#define SAVE_SCALAR_TYPE(POINTER,TYPE,VALUE)    \
  {                                             \
    int n;                                      \
    memcpy(POINTER,&(VALUE),n=sizeof(TYPE));	\
    (POINTER) += n;                             \
  }

#define LOAD_SCALAR_TYPE(POINTER,TYPE,VALUE)    \
  {                                             \
    int n;                                      \
    memcpy(&VALUE,POINTER,n=sizeof(TYPE));	\
    (POINTER) += n;                             \
  }

//--------------------------------------------------

#define SIZE_ARRAY_TYPE(COUNT,TYPE,ARRAY,SIZE)  \
  {                                             \
    (COUNT) += sizeof(int);                     \
    (COUNT) += (SIZE)*sizeof(TYPE);             \
  }

#define SAVE_ARRAY_TYPE(POINTER,TYPE,ARRAY,SIZE)        \
  {                                                     \
    int n,length = (SIZE);                              \
    memcpy(POINTER,&length, n=sizeof(int));             \
    (POINTER) += n;                                     \
    memcpy(POINTER,&ARRAY[0],n=length*sizeof(TYPE));	\
    (POINTER) += n;                                     \
  }
#define LOAD_ARRAY_TYPE(POINTER,TYPE,ARRAY,SIZE)        \
  {                                                     \
    int n,length = (SIZE);                              \
    memcpy(&length, POINTER, n=sizeof(int));		\
    (POINTER) += n;                                     \
    memcpy(&ARRAY[0],POINTER,n=length*sizeof(TYPE));	\
    (POINTER) += n;                                     \
  }

//--------------------------------------------------

#define SIZE_STRING_TYPE(COUNT,STRING)          \
  {                                             \
    (COUNT) += sizeof(int);                     \
    (COUNT) += (STRING).size()*sizeof(char);    \
  }
#define SAVE_STRING_TYPE(POINTER,STRING)                        \
  {                                                             \
    int n,length = (STRING).size();                             \
    memcpy(POINTER,&length, n=sizeof(int));                     \
    (POINTER) += n;                                             \
    memcpy(POINTER,(STRING).data(),n=length*sizeof(char));	\
    (POINTER) += n;                                             \
  }

#define LOAD_STRING_TYPE(POINTER,STRING)                \
  {                                                     \
    int n,length;                                       \
    memcpy(&length, POINTER, n=sizeof(int));            \
    (POINTER) += n;                                     \
    (STRING).resize(length);                            \
    char * string = (char *)(STRING).data();            \
    memcpy(string,POINTER,n=length*sizeof(char));       \
    (POINTER) += n;                                     \
  }

#define SIZE_VECTOR_TYPE(COUNT,TYPE,VECTOR)                     \
  {                                                             \
    (COUNT) += sizeof(int);                                     \
    (COUNT) += sizeof(TYPE)*(VECTOR).size();                    \
  }
#define SAVE_VECTOR_TYPE(POINTER,TYPE,VECTOR)                           \
  {                                                                     \
  int size = (VECTOR).size();                                           \
  memcpy(POINTER,&size, sizeof(int));                                   \
  (POINTER) += sizeof(int);                                             \
  memcpy(POINTER,(TYPE*)&(VECTOR)[0],size*sizeof(TYPE));                \
  (POINTER) += size*sizeof(TYPE);                                       \
  }
#define LOAD_VECTOR_TYPE(POINTER,TYPE,VECTOR)                   \
  {                                                             \
    int size;                                                   \
    memcpy(&size, POINTER, sizeof(int));                        \
    (POINTER) += sizeof(int);                                   \
    (VECTOR).resize(size);                                      \
    memcpy((TYPE*)(VECTOR).data(),POINTER,size*sizeof(TYPE));	\
    (POINTER) += size*sizeof(TYPE);                             \
  }


//--------------------------------------------------

#define SIZE_VECTOR_VECTOR_TYPE(COUNT,TYPE,VECTOR)              \
  {                                                             \
    (COUNT) += sizeof(int);                                     \
    for (std::size_t i=0; i<(VECTOR).size(); i++) {             \
      SIZE_VECTOR_TYPE(COUNT,TYPE,(VECTOR)[i]);                 \
    }                                                           \
  }
#define SAVE_VECTOR_VECTOR_TYPE(POINTER,TYPE,VECTOR)    \
  {                                                     \
    int size = (VECTOR).size();                         \
    memcpy(POINTER,&size, sizeof(int));                 \
    (POINTER) += sizeof(int);                           \
    for (std::size_t i=0; i<(VECTOR).size(); i++) {     \
      SAVE_VECTOR_TYPE(POINTER,TYPE,(VECTOR)[i]);       \
    }                                                   \
  }
#define LOAD_VECTOR_VECTOR_TYPE(POINTER,TYPE,VECTOR)    \
  {                                                     \
    int size;                                           \
    memcpy(&size, POINTER, sizeof(int));                \
    (POINTER) += sizeof(int);                           \
    (VECTOR).resize(size);                              \
    for (std::size_t i=0; i<(VECTOR).size(); i++) {     \
      LOAD_VECTOR_TYPE(POINTER,TYPE,(VECTOR)[i]);       \
    }                                                   \
  }

//--------------------------------------------------

#define SIZE_OBJECT_TYPE(COUNT,OBJECT)          \
  {						\
    (COUNT) += (OBJECT).data_size();            \
  }

#define SAVE_OBJECT_TYPE(POINTER,OBJECT)        \
  {                                             \
    (POINTER) = (OBJECT).save_data(POINTER);    \
  }

#define LOAD_OBJECT_TYPE(POINTER,OBJECT)        \
  {                                             \
    (POINTER) = (OBJECT).load_data(POINTER);    \
  }

//--------------------------------------------------

#define SIZE_OBJECT_PTR_TYPE(COUNT,TYPE,OBJECT_PTR)     \
  {                                                     \
    int have_data = ((OBJECT_PTR) != nullptr);          \
    (COUNT) += sizeof(int);                             \
    if (have_data) {                                    \
      (COUNT) += (OBJECT_PTR)->data_size();             \
    }                                                   \
  }

#define SAVE_OBJECT_PTR_TYPE(POINTER,TYPE,OBJECT_PTR)   \
  {                                                     \
    int have_data = ((OBJECT_PTR) != nullptr);          \
    memcpy ((POINTER),&have_data,sizeof(int));          \
    (POINTER) += sizeof(int);                           \
    if (have_data) {                                    \
      (POINTER) = (OBJECT_PTR)->save_data(POINTER);     \
    }                                                   \
  }

#define LOAD_OBJECT_PTR_TYPE(POINTER,TYPE,OBJECT_PTR)   \
  {                                                     \
    int have_data;                                      \
    memcpy(&have_data,(POINTER),sizeof(int));           \
    (POINTER) += sizeof(int);                           \
    if (have_data) {                                    \
      (OBJECT_PTR) = new TYPE;                          \
      (POINTER) = (OBJECT_PTR)->load_data(POINTER);     \
    } else {                                            \
      (OBJECT_PTR) = nullptr;                           \
    }                                                   \
  }

//--------------------------------------------------

#define SIZE_MAP_TYPE(COUNT,TYPE_1,TYPE_2,MAP)  \
  {                                             \
    (COUNT) += sizeof(int);			\
    (COUNT) += (MAP).size() *                   \
      (sizeof(TYPE_1) + sizeof(TYPE_2));        \
  }
#define SAVE_MAP_TYPE(POINTER,TYPE_1,TYPE_2,MAP)                \
  {                                                             \
    int size = (MAP).size();                                    \
    memcpy(POINTER,&size, sizeof(int));                         \
    (POINTER) += sizeof(int);                                   \
    auto iter = (MAP).begin();                                  \
    while (iter != (MAP).end()) {                               \
      memcpy(POINTER,(TYPE_1*)&iter->first,sizeof(TYPE_1));     \
      (POINTER) += sizeof(TYPE_1);                              \
      memcpy(POINTER,(TYPE_2*)&iter->second,sizeof(TYPE_2));    \
      (POINTER) += sizeof(TYPE_2);                              \
      ++iter;                                                   \
    }                                                           \
  }
#define LOAD_MAP_TYPE(POINTER,TYPE_1,TYPE_2,MAP)        \
  {                                                     \
    int size;                                           \
    memcpy(&size, POINTER, sizeof(int));                \
    (POINTER) += sizeof(int);                           \
    for (int i=0; i<size; i++) {                        \
      TYPE_1 first;                                     \
      TYPE_2 second;                                    \
      memcpy((TYPE_1*)&first,POINTER,sizeof(TYPE_1));   \
      (POINTER) += sizeof(TYPE_1);                      \
      memcpy((TYPE_2*)&second,POINTER,sizeof(TYPE_2));  \
      (POINTER) += sizeof(TYPE_2);                      \
      (MAP)[first] = second;                            \
    }                                                   \
  }

//--------------------------------------------------

#define SIZE_SET_TYPE(COUNT,TYPE,SET)           \
  {                                             \
    (COUNT) += sizeof(int);			\
    (COUNT) += (SET).size() * sizeof(TYPE);     \
  }
#define SAVE_SET_TYPE(POINTER,TYPE,SET)                         \
  {                                                             \
    int size = (SET).size();                                    \
    memcpy(POINTER,&size, sizeof(int));                         \
    (POINTER) += sizeof(int);                                   \
    auto iter = (SET).begin();                                  \
    while (iter != (SET).end()) {                               \
      memcpy(POINTER,(TYPE*)&(*iter),sizeof(TYPE));      \
      (POINTER) += sizeof(TYPE);                                \
      ++iter;                                                   \
    }                                                           \
  }
#define LOAD_SET_TYPE(POINTER,TYPE,SET)                 \
  {                                                     \
    int size;                                           \
    memcpy(&size, POINTER, sizeof(int));                \
    (POINTER) += sizeof(int);                           \
    for (int i=0; i<size; i++) {                        \
      TYPE first;                                       \
      memcpy((TYPE*)&first,POINTER,sizeof(TYPE));       \
      (POINTER) += sizeof(TYPE);                        \
      (SET).insert(first);                              \
    }                                                   \
  }
//--------------------------------------------------

/// Type for CkMyPe(); used for Block() constructor to differentiate
/// from Block(int)
typedef unsigned process_type;

/// Namespace for global constants and functions
namespace cello {

  // Constants

  // pi
  const double pi = 3.14159265358979324;

  // precision functions
  double machine_epsilon     (int);
  template <class T>
  T err_rel (const T & a, const T & b)
  {  return (a != 0.0) ? fabs((a - b) / a) : fabs(a-b);  }

  template <class T>
  T err_abs (const T & a, const T & b)
  {  return fabs(a-b);  }

  template <class T>
  T sum (const T * array,
         int mx, int my, int mz,
         int ox, int oy, int oz,
         int nx, int ny, int nz)
  {
    T s = 0.0;
    int o = ox + mx*(oy + my*oz);
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
        for (int ix=0; ix<nx; ix++) {
          int i=ix+mx*(iy+my*iz) + o;
          s += array[i];
        }
      }
    }
    return s;
  }

  template <class T>
  void copy (T * array_d,
             int mdx, int mdy, int mdz,
             int odx, int ody, int odz,
             const T * array_s,
             int msx, int msy, int msz,
             int osx, int osy, int osz,
             int nx, int ny, int nz)
  {
    int os = osx + msx*(osy + msy*osz);
    int od = odx + mdx*(ody + mdy*odz);
    T sum = 0.0;
    int count=0;
    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
        for (int ix=0; ix<nx; ix++) {
          int is = os + ix+msx*(iy+msy*iz);
          int id = od + ix+mdx*(iy+mdy*iz);
          array_d[id] = array_s[is];
          sum += array_d[id];
          ++count ;
        }
      }
    }
  }

  int digits_max(int precision);

  // type_enum functions (prefered)
  extern bool type_is_float(int type);
  extern bool type_is_int(int type);
  extern bool type_is_valid(int type);

  extern const char * type_name[NUM_TYPES];
  extern const int type_bytes[NUM_TYPES];

  // precision_enum functions (depreciated)

  int sizeof_precision       (precision_type);
  int is_precision_supported (precision_type);
  extern const char * precision_name[7];

  /// converts a precision_enum to type_enum
  ///
  /// at present, this is a trivial operation. This primarily exists to make
  /// the developer's intention clear
  inline int convert_enum_precision_to_type(precision_type precision)
  { return precision; }

  inline void hex_string(char str[], int length)
  {
    //hexadecimal characters
    char hex_characters[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
    static std::default_random_engine generator;
    static std::uniform_int_distribution<int> distribution(0,15);
    int i;
    for(i=0;i<length;i++)
      {
        str[i]=hex_characters[distribution(generator)];
      }
    str[length]=0;
  }

  /// Check a number for zero, NAN, and inf

  template <class T>
  void check (T value, const char * message,
	      const char * file, int line)
  {
    if (std::fpclassify(value) == FP_ZERO) {
      printf ("WARNING: %s:%d %s zero\n", file,line,message);
    }
    if (std::fpclassify(value) == FP_NAN) {
      printf ("WARNING: %s:%d %s nan\n", file,line,message);
    }
    if (std::fpclassify(value) == FP_INFINITE) {
      printf ("WARNING: %s:%d %s inf\n", file,line,message);
    }
#ifdef DEBUG_CHECK
    if (sizeof(value)==sizeof(float))
      CkPrintf ("DEBUG_CHECK %s = %25.15g\n",message,value);
    if (sizeof(value)==sizeof(double))
      CkPrintf ("DEBUG_CHECK %s = %25.15lg\n",message,value);
    if (sizeof(value)==sizeof(long double))
      CkPrintf ("DEBUG_CHECK %s = %25.15Lg\n",message,value);
#endif
  }

  inline int index_static()
  { return CkMyPe() % CONFIG_NODE_SIZE; }

  inline void af_to_xyz (int axis, int face, int r3[3])
  {
    r3[0] = (axis==0) ? 2*face-1 : 0;
    r3[1] = (axis==1) ? 2*face-1 : 0;
    r3[2] = (axis==2) ? 2*face-1 : 0;
  }

  /// Return a pointer to the Simulation object on this process
  Simulation *    simulation();
  /// Return a pointer to the Factory object on this process
  const Factory * factory();
  /// Return a proxy for the Block chare array of Blocks
  CProxy_Block    block_array();
  /// Return a pointor to the Problem object defining the problem being solved
  Problem *       problem();
  /// Return a pointor to the ith Bounday object
  Boundary *      boundary(int i);
  /// Return a pointer to the Hierarchy object defining the mesh hierarchy
  Hierarchy *     hierarchy();
  /// Return a pointer to the Config object containing user parameters values
  const Config *  config();
  /// Return a pointer to the Parameters object
  const Parameters *  parameters();
  /// Return a pointer to the FieldDescr object defining fields on Blocks
  FieldDescr *    field_descr();
  /// Return a pointer to the Groupings object defining field groups
  Grouping *      field_groups();
  /// Define field needed by method or solver
  int define_field (std::string field_name, int cx=0, int cy=0, int cz=0);

  /// Define field needed by method or solver, and also add to a group
  int define_field_in_group (std::string field_name,
                             std::string group_name,
                             int cx=0, int cy=0, int cz=0);
  /// Call after adding all fields, temporary or permanent
  void finalize_fields ();
  /// Return a pointer to the ParticledDescr object defining particles on Blocks
  ParticleDescr * particle_descr();
  /// Return a pointer to the Groupings object defining particle groups
  Grouping *      particle_groups();
  /// Return a pointer to the Monitor object for writing output to stdout
  Monitor *       monitor();
  /// Return a pointer to the Units object
  Units *         units();
  /// Return reference to in indexed Refresh object
  Refresh *       refresh(int ir);
  /// Return the ScalarDescr object defining Block long double Scalar
  /// data values
  ScalarDescr *   scalar_descr_long_double();
  /// Return the ScalarDescr object defining Block double Scalar data values
  ScalarDescr *   scalar_descr_double();
  /// Return the ScalarDescr object defining Block int Scalar data values
  ScalarDescr *   scalar_descr_int();
  /// Return the ScalarDescr object defining Block long long Scalar data values
  ScalarDescr *   scalar_descr_long_long();
  /// Return the ScalarDescr object defining Block Sync counter Scalar
  /// data values
  ScalarDescr *   scalar_descr_sync();
  /// Return the ScalarDescr object defining Block pointer Scalar data
  /// values
  ScalarDescr *   scalar_descr_void();
  /// Return the ScalarDescr object defining Block index Scalar data values
  ScalarDescr *   scalar_descr_index();

  /// Return the ith Output object
  Output *        output (int index);
  /// Return the ith Solver
  Solver *        solver(int index);
  /// Return the dimensional rank of the simulation
  int             rank ();
  /// Return the number of children each Block may have
  int             num_children();
  int             num_children(int rank);
  /// Return the number of Blocks on this process
  size_t          num_blocks_process();
  /// Return the cell volume at the given level relative to the root level
  double          relative_cell_volume (int level);
  //----------------------------------------------------------------------

  /// Return the file name for the format and given arguments
  std::string expand_name
  (const std::vector <std::string> * file_name,
   int counter, int cycle, double time);

  /// Create a directory if it doesn't already exist (returns whether it
  /// already existed)
  bool ensure_directory_exists(const std::string& dir_name);

  /// Return the path to the subdirectory specified by the given format and
  /// arguments this file group output and creates the subdirectory if it
  /// doesn't already exist.
  ///
  /// In the event that path_name doesn't specify a path, the path to the
  /// current directory, ".", is returned
  std::string create_directory
  (const std::vector <std::string> * path_name,
   int counter, int cycle, double time, bool & already_exists);

  //----------------------------------------------------------------------

  // this is a common workaround used to raise a compile-time error in the
  // elsebranch of a constexpr-if statement (as is down directly below
  template<class> inline constexpr bool dummy_false_v_ = false;

  /// returns the type_enum associated with the template type argument ``T``.
  ///
  /// @tparam T The type for which the enum value is returned. Any ``const`` or
  //     ``volatile`` qualifiers are automatically shed from it.
  /// @tparam unknown_on_fail When ``true``, the this returns ``type_unknown``
  ///     when ``T`` isn't recognized. Otherwise, the program fails to compile
  ///     (with a compile-time error), when ``T`` isn't recognized.
  ///
  /// @note
  /// This will not return type_default
  template<typename T, bool unknown_on_fail = false>
  constexpr int get_type_enum() noexcept {
    using T_ = typename std::remove_cv_t<T>;

    if constexpr (std::is_same_v<T_,float>) {
      return type_single;
    } else if constexpr (std::is_same_v<T_, double>) {
      return type_double;
    } else if constexpr (std::is_same_v<T_, long double> && (sizeof(T_)==10)) {
      return type_extended80;
    } else if constexpr (std::is_same_v<T_, long double> && (sizeof(T_)==12)) {
      return type_extended96;
    } else if constexpr (std::is_same_v<T_, long double> && (sizeof(T_)==16)) {
      return type_quadruple;
    } else if constexpr (std::is_same_v<T_, char>) {
      return type_char;
    } else if constexpr (std::is_same_v<T_, short>) {
      return type_short;
    } else if constexpr (std::is_same_v<T_, int>) {
      return type_int;
    } else if constexpr (std::is_same_v<T_, long long>) {
      return type_long_long;
#if 0
    // it's unclear to me at this time if the following is worth including.
    // - On the one hand, fixed-width integer types are not guaranteed to be
    //   defined on all systems.
    // - On the other hand, we could probably provide custom definitions. In
    //   the case where the smallest integer types can't be defined (b/c a
    //   machine's definition of a byte exceeds 8 bits), we will encounter
    //   other problems anyways...
    } else if constexpr (std::is_same_v<T_, std::int8_t>){
      return type_char;
    } else if constexpr (std::is_same_v<T_, std::int16_t>){
      return type_short;
    } else if constexpr (std::is_same_v<T_, std::int32_t>){
      return type_int;
    } else if constexpr (std::is_same_v<T_, std::int64_t>){
      return type_long_long;
#endif
    } else if constexpr (unknown_on_fail) {
      return type_unknown;
    } else {
      static_assert(dummy_false_v_<T_>, "can't convert type!");
    }
  }

  //----------------------------------------------------------------------

  inline int color_get_rgb (std::string color_name) {
    // CSS3 extended standard colors
    // NOTE: defined here instead of e.g. cello::io since accessed in parameters
    // which should not depend on io component
    // source:
    // e.g. https://www.tutorialrepublic.com/css-reference/css-color-names.php
    const static std::map<std::string,int> color_map =
      {
        { "aliceblue" , 0xF0F8FF },
        {  "antiquewhite",   0xFAEBD7 },
        {  "aqua",   0x00FFFF },
        {  "aquamarine",   0x7FFFD4 },
        {  "azure",   0xF0FFFF },
        {  "beige",   0xF5F5DC },
        {  "bisque",   0xFFE4C4 },
        {  "black",   0x000000 },
        {  "blanchedalmond",   0xFFEBCD },
        {  "blue",   0x0000FF },
        {  "blueviolet",   0x8A2BE2 },
        {  "brown",   0xA52A2A },
        {  "burlywood",   0xDEB887 },
        {  "cadetblue",   0x5F9EA0 },
        {  "chartreuse",   0x7FFF00 },
        {  "chocolate",   0xD2691E },
        {  "coral",   0xFF7F50 },
        {  "cornflowerblue",   0x6495ED },
        {  "cornsilk",   0xFFF8DC },
        {  "crimson",   0xDC143C },
        {  "cyan",   0x00FFFF },
        {  "darkblue",   0x00008B },
        {  "darkcyan",   0x008B8B },
        {  "darkgoldenrod",   0xB8860B },
        {  "darkgray",   0xA9A9A9 },
        {  "darkgreen",   0x006400 },
        {  "darkgrey",   0xA9A9A9 },
        {  "darkkhaki",   0xBDB76B },
        {  "darkmagenta",   0x8B008B },
        {  "darkolivegreen",   0x556B2F },
        {  "darkorange",   0xFF8C00 },
        {  "darkorchid",   0x9932CC },
        {  "darkred",   0x8B0000 },
        {  "darksalmon",   0xE9967A },
        {  "darkseagreen",   0x8FBC8F },
        {  "darkslateblue",   0x483D8B },
        {  "darkslategray",   0x2F4F4F },
        {  "darkslategrey",   0x2F4F4F },
        {  "darkturquoise",   0x00CED1 },
        {  "darkviolet",   0x9400D3 },
        {  "deeppink",   0xFF1493 },
        {  "deepskyblue",   0x00BFFF },
        {  "dimgray",   0x696969 },
        {  "dimgrey",   0x696969 },
        {  "dodgerblue",   0x1E90FF },
        {  "firebrick",   0xB22222 },
        {  "floralwhite",   0xFFFAF0 },
        {  "forestgreen",   0x228B22 },
        {  "fuchsia",   0xFF00FF },
        {  "gainsboro",   0xDCDCDC },
        {  "ghostwhite",   0xF8F8FF },
        {  "gold",   0xFFD700 },
        {  "goldenrod",   0xDAA520 },
        {  "gray",   0x808080 },
        {  "green",   0x008000 },
        {  "greenyellow",   0xADFF2F },
        {  "grey",   0x808080 },
        {  "honeydew",   0xF0FFF0 },
        {  "hotpink",   0xFF69B4 },
        {  "indianred",   0xCD5C5C },
        {  "indigo",   0x4B0082 },
        {  "ivory",   0xFFFFF0 },
        {  "khaki",   0xF0E68C },
        {  "lavender",   0xE6E6FA },
        {  "lavenderblush",   0xFFF0F5 },
        {  "lawngreen",   0x7CFC00 },
        {  "lemonchiffon",   0xFFFACD },
        {  "lightblue",   0xADD8E6 },
        {  "lightcoral",   0xF08080 },
        {  "lightcyan",   0xE0FFFF },
        {  "lightgoldenrodyellow",   0xFAFAD2 },
        {  "lightgray",   0xD3D3D3 },
        {  "lightgreen",   0x90EE90 },
        {  "lightgrey",   0xD3D3D3 },
        {  "lightpink",   0xFFB6C1 },
        {  "lightsalmon",   0xFFA07A },
        {  "lightseagreen",   0x20B2AA },
        {  "lightskyblue",   0x87CEFA },
        {  "lightslategray",   0x778899 },
        {  "lightslategrey",   0x778899 },
        {  "lightsteelblue",   0xB0C4DE },
        {  "lightyellow",   0xFFFFE0 },
        {  "lime",   0x00FF00 },
        {  "limegreen",   0x32CD32 },
        {  "linen",   0xFAF0E6 },
        {  "magenta",   0xFF00FF },
        {  "maroon",   0x800000 },
        {  "mediumaquamarine",   0x66CDAA },
        {  "mediumblue",   0x0000CD },
        {  "mediumorchid",   0xBA55D3 },
        {  "mediumpurple",   0x9370DB },
        {  "mediumseagreen",   0x3CB371 },
        {  "mediumslateblue",   0x7B68EE },
        {  "mediumspringgreen",   0x00FA9A },
        {  "mediumturquoise",   0x48D1CC },
        {  "mediumvioletred",   0xC71585 },
        {  "midnightblue",   0x191970 },
        {  "mintcream",   0xF5FFFA },
        {  "mistyrose",   0xFFE4E1 },
        {  "moccasin",   0xFFE4B5 },
        {  "navajowhite",   0xFFDEAD },
        {  "navy",   0x000080 },
        {  "oldlace",   0xFDF5E6 },
        {  "olive",   0x808000 },
        {  "olivedrab",   0x6B8E23 },
        {  "orange",   0xFFA500 },
        {  "orangered",   0xFF4500 },
        {  "orchid",   0xDA70D6 },
        {  "palegoldenrod",   0xEEE8AA },
        {  "palegreen",   0x98FB98 },
        {  "paleturquoise",   0xAFEEEE },
        {  "palevioletred",   0xDB7093 },
        {  "papayawhip",   0xFFEFD5 },
        {  "peachpuff",   0xFFDAB9 },
        {  "peru",   0xCD853F },
        {  "pink",   0xFFC0CB },
        {  "plum",   0xDDA0DD },
        {  "powderblue",   0xB0E0E6 },
        {  "purple",   0x800080 },
        {  "rebeccapurple",   0x663399 },
        {  "red",   0xFF0000 },
        {  "rosybrown",   0xBC8F8F },
        {  "royalblue",   0x4169E1 },
        {  "saddlebrown",   0x8B4513 },
        {  "salmon",   0xFA8072 },
        {  "sandybrown",   0xF4A460 },
        {  "seagreen",   0x2E8B57 },
        {  "seashell",   0xFFF5EE },
        {  "sienna",   0xA0522D },
        {  "silver",   0xC0C0C0 },
        {  "skyblue",   0x87CEEB },
        {  "slateblue",   0x6A5ACD },
        {  "slategray",   0x708090 },
        {  "slategrey",   0x708090 },
        {  "snow",   0xFFFAFA },
        {  "springgreen",   0x00FF7F },
        {  "steelblue",   0x4682B4 },
        {  "tan",   0xD2B48C },
        {  "teal",   0x008080 },
        {  "thistle",   0xD8BFD8 },
        {  "tomato",   0xFF6347 },
        {  "turquoise",   0x40E0D0 },
        {  "violet",   0xEE82EE },
        {  "wheat",   0xF5DEB3 },
        {  "white",   0xFFFFFF },
        {  "whitesmoke",   0xF5F5F5 },
        {  "yellow",   0xFFFF00 },
        {  "yellowgreen",   0x9ACD32 }
      };

    if (color_map.find(color_name) != color_map.end()) {
      return color_map.at(color_name);
    } else {
      return -1;
    }
  }

  /// Returns whether the current cycle is the first cycle after starting the
  /// simulation
  ///
  /// This differentiates between the different kinds of initial cycles
  ///
  /// @note
  /// While the Cello-layer is not actually responsible for implementing
  /// non-charm-based restart functionality, it's useful to store the
  /// functionality in this layer.
  ///
  /// @note
  /// At the time of writing, we can't actually detect when a charm-based
  /// restart has occured - the program aborts with an error in that case
  /// (however, the option exists to make the code more explicit)
  bool is_initial_cycle(InitCycleKind kind) noexcept;
  bool is_initial_cycle(int cycle, InitCycleKind kind) noexcept;
}

#endif /* CELLO_HPP */
