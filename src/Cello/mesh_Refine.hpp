// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Refine.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the Refine Base class
///

#ifndef MESH_REFINE_HPP
#define MESH_REFINE_HPP

class Refine {

  /// @class    Refine
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Refine() throw();

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (FieldBlock * field_block,
		     const FieldDescr * field_descr) throw () = 0;

};

#endif /* MESH_REFINE_HPP */

