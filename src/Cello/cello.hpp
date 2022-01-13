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
#include <vector>

#include <charm++.h>

#include "pup_stl.h"

#include "cello_Sync.hpp"

// #define DEBUG_CHECK

class Block;
class Boundary;
class Config;
class CProxy_Block;
class FieldDescr;
class Field;
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


/// Length of hex message tags used for debugging
#define TAG_LEN 8

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

//--------------------------------------------------

#define SIZE_VECTOR_TYPE(COUNT,TYPE,VECTOR)     \
  {						\
    (COUNT) += sizeof(int);			\
    (COUNT) += sizeof(TYPE)*(VECTOR).size();    \
  }
#define SAVE_VECTOR_TYPE(POINTER,TYPE,VECTOR)                   \
  {                                                             \
    int size = (VECTOR).size();                                 \
    memcpy(POINTER,&size, sizeof(int));                         \
    (POINTER) += sizeof(int);                                   \
    memcpy(POINTER,(TYPE*)&(VECTOR)[0],size*sizeof(TYPE));      \
    (POINTER) += size*sizeof(TYPE);                           \
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

  // ergs per eV
  const double erg_eV = 1.60217653E-12;

  // eV per erg
  const double eV_erg = 6.24150948E11;

  // Boltzman constant in CGS
  const double kboltz = 1.3806504e-16;

  // Solar mass in CGS
  const double mass_solar = 1.98841586e33;

  // Hydrogen mass in CGS
  const double mass_hydrogen = 1.67262171e-24;

  // Electron mass in CGS
  const double mass_electron = 9.10938291E-28;

  // parsec in CGS
  const double pc_cm  = 3.0856775809623245E18;

  // Kiloparsec in CGS
  const double kpc_cm = 3.0856775809623245E21;

  // Megaparsec in CGS
  const double Mpc_cm = 3.0856775809623245E24;

  // speed of light in CGS
  const double clight = 29979245800.0;

  // Gravitational constant in CGS
  const double grav_constant = 6.67384E-8;

  // year in seconds
  const double yr_s = 3.1556952E7;

  // kyr in seconds
  const double kyr_s = 3.1556952E10;

  // Myr in seconds
  const double Myr_s = 3.1556952E13;

  // Approximate mean molecular weight of metals
  const double mu_metal = 16.0;

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
  /// Return the ScalarDescr object defining Block Sync counter Scalar
  /// data values
  ScalarDescr *   scalar_descr_sync();
  /// Return the ScalarDescr object defining Block pointer Scalar data
  /// values
  ScalarDescr *   scalar_descr_void();

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
  (const std::vector <std::string> * file_name, int counter, Block * block);

  /// Return the path for this file group output.  Creates
  /// the subdirectories if they don't exist
  std::string directory
  (const std::vector <std::string> * path_name, int counter, Block * block);

  
}

#endif /* CELLO_HPP */
