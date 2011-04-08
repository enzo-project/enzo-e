// $Id: mesh_PatchMpi.hpp 2181 2011-04-07 00:43:09Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_PATCH_MPI_HPP
#define MESH_PATCH_MPI_HPP

/// @file     mesh_PatchMpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @todo     Move "size" to Block's, since that's Field-centric
/// @brief    [\ref Mesh] Declaration of the interface for the PatchMpi class

#ifndef CONFIG_USE_CHARM

class PatchMpi : public Patch
{

  /// @class    PatchMpi
  /// @ingroup  Mesh
  /// @brief [\ref Mesh] Represent a distributed box of uniform
  /// (non-adaptive) data

public: // interface

  /// Constructor for given Patch size and blocking count
  PatchMpi(Factory * factory,
	GroupProcess * group_process,
	int nx,   int ny,  int nz,
	int nbx,  int nby, int nbz,
	double xm, double ym, double zm,
	double xp, double yp, double zp) throw();

  /// Destructor
  virtual ~PatchMpi() throw();

public: // virtual functions

  /// Allocate local blocks
  virtual void allocate_blocks(FieldDescr * field_descr) throw();

  /// Return the ith local block
  virtual Block * local_block(size_t i) const throw();

  /// Deallocate local blocks
  virtual void deallocate_blocks() throw();

protected: // attributes

  /// Array of blocks ib associated with this process
  std::vector<Block * > block_;
};

#endif /* ! CONFIG_USE_CHARM */

#endif /* MESH_PATCH_MPI_HPP */

