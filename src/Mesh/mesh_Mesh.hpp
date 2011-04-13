// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Mesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 10 15:38:40 PST 2009
/// @brief    [\ref Mesh] Declaration of the Mesh class

#ifndef MESH_MESH_HPP
#define MESH_MESH_HPP

class Factory;

class Mesh {

  /// @class    Mesh
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Adaptive mesh refinement hierarchy

public: // interface

  /// Initialize an Mesh object
  Mesh(const Factory * factory,
       int nx,  int ny,  int nz,
       int nbx, int nby, int nbz) throw ();

  /// Delete the Mesh object
  virtual ~Mesh() throw ();

  //----------------------------------------------------------------------

  /// Set domain lower extent
  void set_lower(double nx, double ny, double nz) throw ();

  /// Set domain upper extent
  void set_upper(double nx, double ny, double nz) throw ();
  
  //----------------------------------------------------------------------

  /// Return dimension
  int dimension() const throw ();

  /// Return domain lower extent
  void lower(double * nx, double * ny = 0, double * nz = 0) const throw ();

  /// Return domain upper extent
  void upper(double * nx, double * ny = 0, double * nz = 0) const throw ();

  //----------------------------------------------------------------------

  /// Return the total number of local patches
  size_t num_patches() const throw();

  /// Return the ith patch.  Patch[0] is the root
  Patch * patch(size_t i) throw();

  /// Return the ith patch
  Patch * patch(size_t i) const throw();

  /// Insert the given Patch into the list of patches
  virtual void insert_patch(Patch *) throw();

  /// Return the factory object associated with the Mesh
  const Factory * factory () const throw()
  { return factory_; }

protected: // attributes

  /// Factory for creating Simulations, Meshes, Patches and Blocks
  /// [abstract factory design pattern]
  const Factory * factory_;

  /// List of local patchs

  std::vector<Patch *> patch_list_;

  /// Tree defining the MESH hierarchy topology
  //  strict_auto_ptr<TreeK> tree_;
  TreeK * tree_;

  /// Lower extent of the patch
  double lower_[3];

  /// Upper extent of the patch
  double upper_[3];

};

  // /// Set max_level
  // void set_max_level(int max_level) throw ();

  // /// Set refinement factor
  // void set_refine_factor(int refine) throw ();

  // /// Set whether to avoid level jumps
  // void set_balanced(bool balanced) throw ();

  // /// Set whether to backfill levels
  // void set_backfill(bool backfill) throw ();

  // /// Set whether to coalesce patches
  // void set_coalesce(bool coalesce) throw ();
  // //--------------------------------------------------
  // /// Return max_level
  // int max_level() const throw ();

  // /// Return refinement factor
  // int refine_factor() const throw ();

  // /// Return whether to avoid level jumps
  // bool balanced() const throw ();

  // /// Return whether to backfill levels
  // bool backfill() const throw ();

  // /// Return whether to coalesce patches
  // bool coalesce() const throw ();

  // private:
  // /// Refinement factor = 2, 4, etc.
  // /// Parameter Mesh::refine
  // int refine_;

  // /// Maximum level for the hierarchy (0 = unigrid) assuming r=2
  // /// Parameter Mesh::max_level
  // int max_level_;

  // /// Whether the tree is balanced or "full"
  // /// Parameter Mesh::balanced
  // bool balanced_;

  // /// Whether to backfill for refine > 2 to regain r == 2 balance
  // /// Parameter Mesh::backfill
  // bool backfill_;

  // /// Whether to coalesce small patches into one big one
  // /// Parameter Mesh::coalesce
  // bool coalesce_;

#endif /* MESH_MESH_HPP */

