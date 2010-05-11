// $Id: mesh_mesh.hpp 1259 2010-03-02 03:12:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_mesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 10 15:38:40 PST 2009
/// @brief    Declaration of the Mesh class

#ifndef MESH_MESH_HPP
#define MESH_MESH_HPP

#include <memory>
#include <vector>
#include "strict_auto_ptr.h"

class Mesh {

  /// @class    Mesh
  /// @ingroup  Mesh
  /// @brief    Adaptive mesh refinement hierarchy

public: // interface

  /// Initialize an Mesh object
  Mesh() :
    tree_(0),
    dimension_(0),
    max_level_(0),
    refine_(2),
    root_size_(),
    min_patch_size_(4),
    max_patch_size_(-1),
    balanced_(true),
    backfill_(true),
    coalesce_(true)
  {};

  /// Return dimension
  int dimension() 
  { return dimension_; };

  /// Set dimension
  void set_dimension(int dimension) 
  {dimension_ = dimension; };

  /// Return max_level
  int max_level() 
  { return max_level_; };

  /// Set max_level
  void set_max_level(int max_level) 
  {max_level_ = max_level; };

  /// Return max_level
  int refine() 
  {return refine_; };

  void set_refine(int refine) 
  {refine_ = refine; }; 

  /// Return root_size
  std::vector<int> root_size()
  { return root_size_; };

  void set_root_size(std::vector<int> root_size) 
  { root_size_ = root_size; };

  /// Return min_patch_size
  int min_patch_size()
  {return min_patch_size_; };

  void set_min_patch_size(int min_patch_size)
  { min_patch_size_ = min_patch_size; };

  /// Return max_patch_size
  int max_patch_size()
  {return max_patch_size_; };

  void set_max_patch_size(int max_patch_size)
  { max_patch_size_ = max_patch_size; };

  /// Return balanced
  bool balanced()
  {return balanced_; };

  void set_balanced(bool balanced)
  { balanced_ = balanced; };

  /// Return backfill
  bool backfill()
  {return backfill_; };

  void set_backfill(bool backfill)
  { backfill_ = backfill; };

  /// Return coalesce
  bool coalesce()
  {return coalesce_; };

  void set_coalesce(bool coalesce)
  { coalesce_ = coalesce; };
  
private: // attributes

  /// Tree defining the MESH hierarchy topology
  strict_auto_ptr<TreeK> tree_;

  /// Spacial dimensions of the Mesh: 1, 2, or 3
  /// Determined from Domain::extent
  int dimension_;

  /// Maximum level for the hierarchy (0 = unigrid) assuming r=2
  /// Parameter Mesh::max_level
  int max_level_;

  /// Refinement factor = 2, 4, etc.
  /// Parameter Mesh::refine
  int refine_;

  /// Size of root-level mesh.  Move to Field?
  /// Parameter Mesh::root_size
  std::vector<int> root_size_; 

  /// Minimum allowed patch size (Move to Field?)
  /// Parameter Mesh::min_patch_size
  int min_patch_size_;

  /// Maximum allowed patch size (Move to Field?)
  /// Parameter Mesh::max_patch_size
  int max_patch_size_;

  /// Whether the tree is balanced or "full"
  /// Parameter Mesh::balanced
  bool balanced_;

  /// Whether to backfill for refine > 2 to regain r == 2 balance
  /// Parameter Mesh::backfill
  bool backfill_;

  /// Whether to coalesce small patches into one big one
  /// Parameter Mesh::coalesce
  bool coalesce_;

};

#endif /* MESH_MESH_HPP */

