// $Id: mesh_mesh.hpp 1259 2010-03-02 03:12:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_mesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 10 15:38:40 PST 2009
/// @brief    Declaration of the Mesh class

#ifndef MESH_MESH_HPP
#define MESH_MESH_HPP

#include <memory>
#include "strict_auto_ptr.h"

class Mesh {

  /// @class    Mesh
  /// @ingroup  Mesh
  /// @brief    Adaptive mesh refinement hierarchy

public: // interface

  /// Initialize an Mesh object
  Mesh() :
    max_level_(0),
    tree_(0)
  {};

private: // attributes

  /// Maximum level for the hierarchy (0 = unigrid)
  int max_level_;

  /// Tree defining the MESH hierarchy topology
  strict_auto_ptr<TreeK> tree_;


};

#endif /* MESH_MESH_HPP */

