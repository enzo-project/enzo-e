// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the RefineSlope class
///

#ifndef MESH_REFINE_SLOPE_HPP
#define MESH_REFINE_SLOPE_HPP

class RefineSlope {

  /// @class    RefineSlope
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  RefineSlope() throw();

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (FieldBlock * field_block,
		     const FieldDescr * field_descr) throw();

};

#endif /* MESH_REFINE_SLOPE_HPP */

