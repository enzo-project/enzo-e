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

  friend class IoHierarchy;

public: // interface

  /// Initialize a Hierarchy object
  Hierarchy ( const Factory * factory,
	      int dimension, int refinement) throw ();

  /// Delete the Hierarchy object
  virtual ~Hierarchy() throw ();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change

    p | *factory_;
    p | patch_count_;
    p | *patch_tree_;
    PUParray(p,root_size_,3);
    PUParray(p,lower_,3);
    PUParray(p,upper_,3);

  }
#endif

  //----------------------------------------------------------------------

  /// Set domain lower extent
  void set_lower(double x, double y, double z) throw ();

  /// Set domain upper extent
  void set_upper(double x, double y, double z) throw ();
  
  /// Set root-level grid size
  void set_root_size(int nx, int ny, int nz) throw ();

  //----------------------------------------------------------------------

  /// Return dimension
  int dimension() const throw ();

  /// Return domain lower extent
  void lower(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return domain upper extent
  void upper(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return root-level grid size
  void root_size(int * nx, int * ny = 0, int * nz = 0) const throw ();

  //----------------------------------------------------------------------

  /// Return the total number of local patches
  size_t num_patches() const throw();

  /// Return the ith patch.  Patch[0] is the root
  Patch * patch(size_t i) throw();

  /// Return the ith patch
  Patch * patch(size_t i) const throw();

  /// Create the initial root patch
  void create_root_patch (FieldDescr   * field_descr,
			  int nx, int ny, int nz,
			  int nbx, int nby, int nbz,
			  bool allocate_blocks  = true,
			  int process_first     = 0, 
			  int process_last_plus = -1) throw();

  /// Return the factory object associated with the Hierarchy
  const Factory * factory () const throw()
  { return factory_; }

protected: // attributes

  /// Factory for creating Simulations, Hierarchies, Patches and Blocks
  /// [abstract factory design pattern]
  Factory * factory_;

  /// Number of patches (redundant with patch_tree_)
  int patch_count_;

  /// List of local patches
  Tree * patch_tree_;

  /// Size of the root grid
  int root_size_[3];

  /// Lower extent of the hierarchy
  double lower_[3];

  /// Upper extent of the hierarchy
  double upper_[3];

};


#endif /* MESH_HIERARCHY_HPP */

