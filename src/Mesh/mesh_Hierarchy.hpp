// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Hierarchy.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Nov 10 15:38:40 PST 2009
/// @brief    [\ref Mesh] Declaration of the Hierarchy class

#ifndef MESH_HIERARCHY_HPP
#define MESH_HIERARCHY_HPP

class Factory;

class Hierarchy {

  /// @class    Hierarchy
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Adaptive mesh refinement hierarchy

public: // interface

  /// Initialize a Hierarchy object
  Hierarchy(const Factory * factory) throw ();

  /// Delete the Hierarchy object
  virtual ~Hierarchy() throw ();

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

  /// Create the initial root patch
  void create_root_patch (GroupProcess * group_process,
			  FieldDescr   * field_descr,
			  int nx, int ny, int nz,
			  int nbx, int nby, int nbz) throw();

  /// Return the factory object associated with the Hierarchy
  const Factory * factory () const throw()
  { return factory_; }

protected: // attributes

  /// Factory for creating Simulations, Hierarchies, Patches and Blocks
  /// [abstract factory design pattern]
  const Factory * factory_;

  /// List of local patchs

  std::vector<Patch *> patch_list_;

  /// Tree defining the mesh hierarchy topology
  //  strict_auto_ptr<TreeK> tree_;
  TreeK * tree_;

  /// Lower extent of the hierarchy
  double lower_[3];

  /// Upper extent of the hierarchy
  double upper_[3];

};

#endif /* MESH_HIERARCHY_HPP */

