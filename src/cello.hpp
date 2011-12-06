// See LICENSE_CELLO file for license and copyright information

#ifndef CELLO_HPP
#define CELLO_HPP

/// @file    cello.hpp
/// @author  James Bordner (jobordner@ucsd.edu)
/// @date    Thu Nov 11 17:08:38 PST 2010
/// @todo    Move TEMP_CLEAR_VALUE to Memory parameter
/// @brief   Include Cello global configuration settings

//----------------------------------------------------------------------
// COMMON FUNCTIONS
//----------------------------------------------------------------------

template <class T>
inline T MIN(const T &a, const T &b) 
{  return a < b ? a : b; }

template <class T>
inline T MAX(const T &a, const T &b) 
{  return a > b ? a : b; }

inline int INDEX(int ix,int iy,int iz,int nx,int ny) 
{  return ix+nx*(iy+ny*iz); }

//----------------------------------------------------------------------
// GLOBAL DEFINES
//----------------------------------------------------------------------

// Value used to initialize new fields for debugging

#define TEMP_CLEAR_VALUE (-std::numeric_limits<float>::max())

// (std::numeric_limits<float>::max())

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

/// @enum     reduce_enum
/// @brief    Reduction operator, used for image projections
enum reduce_enum {
  reduce_unknown, /// Unknown reduction
  reduce_min,     /// Minimal value along the axis
  reduce_max,     /// Maximal value along the axis
  reduce_avg,     /// Average value along the axis
  reduce_sum      /// Sum of values along the axis
};

//======================================================================
// COMPONENTS
//======================================================================

enum component_enum {
  component_undefined,
  component_disk,
  component_first = component_disk,
  component_error,
  component_field,
  component_memory,
  component_mesh,
  component_method,
  component_monitor,
  component_parallel,
  component_parameters,
  component_particles,
  component_performance,
  component_portal,
  component_simulation,
  num_components = component_simulation
};

extern const char * component_name [];

//======================================================================
// PRECISION
//======================================================================

enum precision_enum {
  // @@@ KEEP IN SYNCH WITH precision_name in cello_precision.cpp
  precision_unknown,     //  unknown precision
  precision_default,     //  default precision
  precision_single,      //  32-bit field data
  precision_double,      //  64-bit field data
  precision_extended80,  //  80-bit field data
  precision_extended96,  //  96-bit field data
  precision_quadruple,   // 128-bit field data
};

#ifdef CONFIG_PRECISION_SINGLE

#   define default_precision precision_single
#   define scalar_type_enzo_float scalar_type_float

#   define SCALAR_DEFINED

#endif

#ifdef CONFIG_PRECISION_DOUBLE

#   define default_precision precision_double
#   define scalar_type_enzo_float scalar_type_double

#   ifdef SCALAR_DEFINED
#      error Both CONFIG_PRECISION_SINGLE and CONFIG_PRECISION_DOUBLE defined
#   endif

#   define SCALAR_DEFINED

#endif

#ifndef SCALAR_DEFINED

#   error Neither CONFIG_PRECISION_SINGLE nor CONFIG_PRECISION_DOUBLE defined

#endif

namespace cello {

  int sizeof_precision       (enum precision_enum);
  int is_precision_supported (enum precision_enum);
  double machine_epsilon     (enum precision_enum);
  extern const char * precision_name[7];

  template <class T>
  T err_rel (T a, T b)
  {  return (a != 0.0) ? fabs((a - b) / a) : fabs(a-b);  }

  template <class T>
  T err_abs (T a, T b)
  {  return fabs(a-b);  };

};

#endif /* CELLO_HPP */
