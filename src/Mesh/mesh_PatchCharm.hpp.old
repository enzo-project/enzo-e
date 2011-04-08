// $Id: mesh_PatchCharm.hpp 2181 2011-04-07 00:43:09Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_PATCH_CHARM_HPP
#define MESH_PATCH_CHARM_HPP

/// @file     mesh_PatchCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @todo     Move "size" to Block's, since that's Field-centric
/// @brief    [\ref Mesh] Declaration of the interface for the PatchCharm class

#ifdef CONFIG_USE_CHARM

#include "enzo.decl.h"

class PatchCharm : public Patch
{

  /// @class    PatchCharm
  /// @ingroup  Mesh
  /// @brief [\ref Mesh] Represent a distributed box of uniform
  /// (non-adaptive) data

public: // interface

  /// Constructor for given PatchCharm size and blocking count
  PatchCharm(Factory * factory,
	GroupProcess * group_process,
	int nx,   int ny,  int nz,
	int nbx,  int nby, int nbz,
	double xm, double ym, double zm,
	double xp, double yp, double zp) throw();

  virtual ~PatchCharm() throw();

public: // virtual functions

  /// Allocate local blocks
  virtual void allocate_blocks(FieldDescr * field_descr) throw();

  /// Return the ith local block
  virtual Block * local_block(size_t i) const throw();

  /// Deallocate local blocks
  virtual void deallocate_blocks() throw();

protected: // attributes

  CProxy_BlockCharm block_;
};

#endif /* CONFIG_USE_CHARM */

#endif /* MESH_PATCH_CHARM_HPP */

