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

  /// Initialize an Hierarchy object
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


  /*
    @@@@@@@@@@@@@@@@@@@@@@@@@@
  /// Set max_level
  void set_max_level(int max_level) throw ();

  /// Return max_level
  int max_level() const throw ();

  /// Set refinement factor
  void set_refine_factor(int refine) throw ();

  /// Return refinement factor
  int refine_factor() const throw ();

  /// Set whether to avoid level jumps
  void set_balanced(bool balanced) throw ();

  /// Return whether to avoid level jumps
  bool balanced() const throw ();

  /// Set whether to backfill levels
  void set_backfill(bool backfill) throw ();

  /// Return whether to backfill levels
  bool backfill() const throw ();

  /// Set whether to coalesce patches
  void set_coalesce(bool coalesce) throw ();

  /// Return whether to coalesce patches
  bool coalesce() const throw ();

    @@@@@@@@@@@@@@@@@@@@@@@@@@
  */

  //----------------------------------------------------------------------
  // I/O
  //----------------------------------------------------------------------

  /// Open a file for the Hierarchy
  void open (File * file, const char * filename, const char * mode) const throw();

  /// Close a file for the Hierarchy
  void close (File * file) const throw();

  /// Read "metadata" or field data associated with the Hierarchy
  void read (File * file, file_content_type file_content) throw ();

  /// Write "metadata" or field data associated with the Hierarchy
  void write(File * file, file_content_type file_content) const throw ();

protected: // attributes

  /// Factory for creating Simulations, Hierarchies, Patches and Blocks
  /// [abstract factory design pattern]
  const Factory * factory_;

  /// List of local patchs

  std::vector<Patch *> patch_list_;

  /// Tree defining the mesh hierarchy topology
  //  strict_auto_ptr<TreeK> tree_;
  TreeK * tree_;

  /// Lower extent of the patch
  double lower_[3];

  /// Upper extent of the patch
  double upper_[3];

  /*
  /// Refinement factor = 2, 4, etc.
  /// Parameter Hierarchy::refine
  int refine_;

  /// Maximum level for the hierarchy (0 = unigrid) assuming r=2
  /// Parameter Hierarchy::max_level
  int max_level_;

  /// Whether the tree is balanced or "full"
  /// Parameter Hierarchy::balanced
  bool balanced_;

  /// Whether to backfill for refine > 2 to regain r == 2 balance
  /// Parameter Hierarchy::backfill
  bool backfill_;

  /// Whether to coalesce small patches into one big one
  /// Parameter Hierarchy::coalesce
  bool coalesce_;

    @@@@@@@@@@@@@@@@@@@@@@@@@@
  */

};

#endif /* MESH_HIERARCHY_HPP */

