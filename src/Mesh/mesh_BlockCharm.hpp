// $Id: mesh_BlockCharm.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_BLOCK_CHARM_HPP
#define MESH_BLOCK_CHARM_HPP

/// @file     mesh_BlockCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:01:51 PST 2011
/// @todo     Change public attributes to private
/// @todo     Dynamically allocate arrays
/// @brief    [\ref ] Declaration of the BlockCharm class

#ifdef CONFIG_USE_CHARM

#include "enzo.decl.h"

class BlockCharm : public Block , public CBase_BlockCharm {

  /// @class    BlockCharm
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] A BlockCharm is a Block for use with CHARM++

public: // interface

  /// Initialize the BlockCharm chare
  BlockCharm
  ( int nx, int ny, int nz,
    double xm, double ym, double zm,
    double hx, double hy, double hz,
    int num_field_blocks);

  /// Migration constructor (e.g. p.27 The CHARM Programming Language
  /// Manual v6.0)
  BlockCharm (CkMigrateMessage *m) {TRACE("Oops")};

  BlockCharm () {TRACE("Oops")};

  virtual ~BlockCharm() throw() {TRACE("Oops")};

public: // charm entry functions

  /// Initialize block for the simulation.
  void p_initial();

  /// Cycle output (Disk and Monitor), compute timestep and stopping
  // /// criteria, exit simulation if done
  // void p_next();

  // /// Update boundary conditions and refresh ghost zones
  // void p_prepare();

  // /// Apply methods to blocks
  // void p_compute();

};

#endif /* CONFIG_USE_CHARM */

#endif /* MESH_BLOCK_CHARM_HPP */

