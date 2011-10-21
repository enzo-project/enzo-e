// See LICENSE_CELLO file for license and copyright information

#ifndef CELLO_HPP
#define CELLO_HPP

/// @file    cello.hpp
/// @author  James Bordner (jobordner@ucsd.edu)
/// @date    Thu Nov 11 17:08:38 PST 2010
/// @todo    Need face_axis_enum?
/// @brief   Include Cello global configuration settings

#include "cello_macros.hpp"
#include "cello_precision.hpp"

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

#define TEMP_CLEAR_VALUE (std::numeric_limits<float>::min())



/// @enum     face_axis_enum
/// @brief    Face [lower|upper][x|y|z]
enum face_axis_enum {
  face_lower_axis_x = 0,
  face_upper_axis_x = 1,
  face_lower_axis_y = 2,
  face_upper_axis_y = 3,
  face_lower_axis_z = 4,
  face_upper_axis_z = 5,
  face_axis_all
};

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

//----------------------------------------------------------------------
/// @enum     reduce_enum
/// @brief    Reduction operator, used for image projections

enum reduce_enum {
  reduce_unknown, /// Unknown reduction
  reduce_min,     /// Minimal value along the axis
  reduce_max,     /// Maximal value along the axis
  reduce_avg,     /// Average value along the axis
  reduce_sum      /// Sum of values along the axis
};

/*********************************************************************
 * COMPONENTS
 **********************************************************************/

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

#endif /* CELLO_HPP */
