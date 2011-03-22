// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_PATCH_HPP
#define MESH_PATCH_HPP

/// @file     mesh_Patch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @todo     Move "size" to Block's, since that's Field-centric
/// @brief    [\ref Mesh] Declaration of the interface for the Patch class

class Mesh;

#include PARALLEL_CHARM_INCLUDE(Patch.decl.h)

PARALLEL_CLASS_DECL(Patch)
{

  /// @class    Patch
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Represent a distributed box of uniform (non-adaptive) data

 public: // interface

#ifdef CONFIG_USE_CHARM
  /// Constructor for CHARM++
  Patch(int nx,   int ny,  int nz,
	int nbx,  int nby, int nbz);
#endif

  /// Constructor for given Patch size and blocking count
  Patch(Factory * factory,
	GroupProcess * group_process,
	int nx,   int ny,  int nz,
	int nbx,  int nby, int nbz) throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  virtual ~Patch() throw();

  /// Copy constructor
  Patch(const Patch & patch,
	FieldDescr * field_descr) throw();

  // /// Assignment operator
  // Patch & operator= (const Patch & patch) throw();

  //----------------------------------------------------------------------

  /// Set domain lower extent
  void set_lower(double x, double y, double z) throw ();

  /// Set domain upper extent
  void set_upper(double x, double y, double z) throw ();
  
  //----------------------------------------------------------------------

  /// Return the size of the patch in number of grid cells
  void size (int * nx, int * ny=0, int * nz=0) const throw();

  /// Return the number of blocks along each dimension
  void blocking (int * nbx, int * nby=0, int * nbz=0) const throw();

  /// Return the layout of the patch, describing processes and blocking
  Layout * layout () const throw();

  /// Return domain lower extent
  void lower(double * x, double * y=0, double * z=0) const throw ();

  /// Return domain upper extent
  void upper(double * x, double * y=0, double * z=0) const throw ();

  //--------------------------------------------------

  /// Allocate local blocks
  void allocate_blocks(FieldDescr * field_descr) throw();

  /// Deallocate local blocks
  void deallocate_blocks() throw();

  /// Whether local blocks are allocated
  bool blocks_allocated() const throw();

  /// Return the total number of local blocks
  size_t num_blocks() const throw();

  /// Return the ith local block
  Block * local_block(size_t i) const throw();

  GroupProcess * group()  const throw()
  { return group_process_; };

 public: // entry functions

#ifdef CONFIG_USE_CHARM

#endif

  //--------------------------------------------------

 protected: // attributes

  /// Factory object for creating Blocks
  Factory * factory_;

  /// Parallel Group for distributing the Mesh across processors
  GroupProcess * group_process_;

  /// Layout: describes blocking, processor range, and block-processor mapping 
  Layout * layout_;

  /// Array of blocks ib associated with this process
  std::vector<Block * > block_;

  /// Size of the patch
  int size_[3];

  /// How the Patch is distributed into Blocks
  int blocking_[3];

  /// Lower extent of the patch
  double lower_[3];

  /// Upper extent of the patch
  double upper_[3];


};

#endif /* MESH_PATCH_HPP */

