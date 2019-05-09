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

#include <stdio.h>
#include <math.h>

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <execinfo.h>

#include <charm++.h>

#include "pup_stl.h"

// #define DEBUG_CHECK

class Config;
class CProxy_Block;
class FieldDescr;
class Hierarchy;
class Monitor;
class Output;
class ParticleDescr;
class Problem;
class ScalarDescr;
class Simulation;
class Solver;
class Units;

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
#endif
#ifdef CONFIG_PRECISION_DOUBLE
#   ifdef SCALAR_DEFINED
#      define SCALAR_ERROR
#   endif
#   define default_precision precision_double
#   define default_type      type_double
#   define SCALAR_DEFINED
#endif
#ifdef CONFIG_PRECISION_QUAD
#   ifdef SCALAR_DEFINED
#      define SCALAR_ERROR
#   endif
#   define default_precision precision_quad
#   define default_type      type_quad
#   define SCALAR_DEFINED
#endif

#ifndef SCALAR_DEFINED
#   error None of CONFIG_PRECISION_[SINGLE|DOUBLE|QUAD] defined
#endif

#ifdef ERROR_SCALAR
#   error Multiple CONFIG_PRECISION_[SINGLE|DOUBLE|QUAD] defined
#endif

/// Macros for sizing, saving, and restoring data from buffers

#define SIZE_ARRAY(COUNT,LIST)			\
  {						\
    (*COUNT) += sizeof(int);			\
    (*COUNT) += sizeof(int)*LIST.size();	\
  }
#define SIZE_VALUE(COUNT,VALUE)			\
  {						\
    (*COUNT) += sizeof(int);			\
  }

#define SAVE_ARRAY(PTR,LIST)				\
  {							\
    int length = LIST.size();				\
    int n;						\
    memcpy((*PTR),&length, n=sizeof(int));		\
    (*PTR)+=n;						\
    memcpy((*PTR),&LIST[0],n=length*sizeof(int));	\
    (*PTR)+=n;						\
  }
#define SAVE_VALUE(PTR,VALUE)			\
  {						\
    int n;					\
    memcpy((*PTR),&VALUE,n=sizeof(int));	\
    (*PTR)+=n;					\
  }

#define LOAD_ARRAY(PTR,LIST)				\
  {							\
    int length;						\
    int n;						\
    memcpy(&length, (*PTR), n=sizeof(int));		\
    (*PTR)+=n;						\
    LIST.resize(length);				\
    memcpy(&LIST[0],(*PTR),n=length*sizeof(int));	\
    (*PTR)+=n;						\
  }
#define LOAD_VALUE(PTR,VALUE)			\
  {						\
    int n;					\
    memcpy(&VALUE,(*PTR),n=sizeof(int));	\
    (*PTR)+=n;					\
  }

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

  // type_enum functions (prefered)
  int sizeof_type (int);
  int is_type_supported (int);

  extern bool type_is_float(int type);
  extern bool type_is_int(int type);
  extern bool type_is_valid(int type);

  extern const char * type_name[NUM_TYPES];
  extern const int type_bytes[NUM_TYPES];

  // precision_enum functions (depreciated)

  int sizeof_precision       (precision_type);
  int is_precision_supported (precision_type);
  extern const char * precision_name[7];


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

  void backtrace(const char * msg);

  inline int index_static()
  { return CkMyPe() % CONFIG_NODE_SIZE; }

  /// Return a pointer to the Simulation object on this process
  Simulation *    simulation();
  /// Return a proxy for the Block chare array of Blocks
  CProxy_Block    block_array();
  /// Return a pointor to the Problem object defining the problem being solved
  Problem *       problem();
  /// Return a pointer to the Hierarchy object defining the mesh hierarchy
  Hierarchy *     hierarchy();
  /// Return a pointer to the Config object containing user parameters values
  const Config *  config();
  /// Return a pointer to the FieldDescr object defining fields on Blocks
  FieldDescr *    field_descr();
  /// Return a pointer to the ParticledDescr object defining particles on Blocks
  ParticleDescr * particle_descr();
  /// Return a pointer to the Monitor object for writing output to stdout
  Monitor *       monitor();
  /// Return a pointer to the Units object
  Units *         units();

  /// Return the ScalarDescr object defining Block long double Scalar data values
  ScalarDescr *   scalar_descr_long_double();
  /// Return the ScalarDescr object defining Block double Scalar data values
  ScalarDescr *   scalar_descr_double();
  /// Return the ScalarDescr object defining Block int Scalar data values
  ScalarDescr *   scalar_descr_int();
  /// Return the ScalarDescr object defining Block Sync counter Scalar data values
  ScalarDescr *   scalar_descr_sync();
  /// Return the ScalarDescr object defining Block pointer Scalar data values
  ScalarDescr *   scalar_descr_void();

  /// Return the ith Output object
  Output *        output (int index);
  /// Return the ith Solver
  Solver *        solver(int index);
  /// Return the dimensional rank of the simulation
  int             rank ();
  /// Return the number of children each Block may have
  int             num_children();
}

#endif /* CELLO_HPP */
