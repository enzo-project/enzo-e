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
#include "pup_stl.h"

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
typedef int face_type;

/// @enum     axis_enum
/// @brief    Axis [x|y|z]
enum axis_enum {
  axis_x = 0,
  axis_y = 1,
  axis_z = 2,
  axis_all
};
typedef int axis_type;

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

// extern const char * component_name [];
/// @enum precision_enum
/// @brief list of known floating-point precision
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
#   define default_precision precision_single
#   define SCALAR_DEFINED
#endif
#ifdef CONFIG_PRECISION_DOUBLE
#   ifdef SCALAR_DEFINED
#      define SCALAR_ERROR
#   endif
#   define default_precision precision_double
#   define SCALAR_DEFINED
#endif
#ifdef CONFIG_PRECISION_QUAD
#   ifdef SCALAR_DEFINED
#      define SCALAR_ERROR
#   endif
#   define default_precision precision_quad
#   define SCALAR_DEFINED
#endif

#ifndef SCALAR_DEFINED
#   error Neither CONFIG_PRECISION_SINGLE nor CONFIG_PRECISION_DOUBLE defined
#endif

#ifdef ERROR_SCALAR
#   error Multiple CONFIG_PRECISION_<SIZE> defined
#endif

/// Namespace for global constants and functions
namespace cello {
  
  // Constants

  const double pi = 3.14159265358979324;
  const double G_cgs  = 6.67384e-8;

  // Functions
  int sizeof_precision       (precision_type);
  int is_precision_supported (precision_type);
  double machine_epsilon     (precision_type);
  extern const char * precision_name[7];

  template <class T>
  T err_rel (const T & a, const T & b)
  {  return (a != 0.0) ? fabs((a - b) / a) : fabs(a-b);  }

  template <class T>
  T err_abs (const T & a, const T & b)
  {  return fabs(a-b);  };

};

#endif /* CELLO_HPP */
