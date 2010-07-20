// $Id: mesh_Patch.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_PATCH_HPP
#define MESH_PATCH_HPP

/// @file     mesh_Patch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the interface for the Patch class

#include "data.hpp"
#include "parallel.hpp"

class Patch {

  /// @class    Patch
  /// @ingroup  Mesh
  /// @brief    Brief description of class Patch.

public: // interface

  /// Constructor
  Patch() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Patch() throw();

  /// Copy constructor
  Patch(const Patch & patch) throw();

  /// Assignment operator
  Patch & operator= (const Patch & patch) throw();

private: // attributes

  /// Parallel layout descriptor
  ParallelLayout * layout_;

  /// Array of local data blocks associated with the patch
  DataBlock * data_block_;

};

#endif /* MESH_PATCH_HPP */

