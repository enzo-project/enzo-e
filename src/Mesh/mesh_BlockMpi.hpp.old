// $Id: mesh_BlockMpi.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_BLOCK_MPI_HPP
#define MESH_BLOCK_MPI_HPP

/// @file     mesh_BlockMpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:01:51 PST 2011
/// @todo     Change public attributes to private
/// @todo     Dynamically allocate arrays
/// @brief    [\ref ] Declaration of the BlockMpi class

#ifndef CONFIG_USE_CHARM

class BlockMpi : public Block {

  /// @class    BlockMpi
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] A BlockMpi is a Block for use with MPI++

public: // interface

  BlockMpi(int ix, int iy, int iz,
	    int nx, int ny, int nz,
	    double xm, double ym, double zm,
	    double xp, double yp, double zp,
	    int num_field_blocks) throw();

  virtual ~BlockMpi() throw (){};

public: // mpi entry functions

  /// Cycle output (Disk and Monitor), compute timestep and stopping
  // /// criteria, exit simulation if done
  // void p_next();

  // /// Update boundary conditions and refresh ghost zones
  // void p_prepare();

  // /// Apply methods to blocks
  // void p_compute();

};

#endif /* ! CONFIG_USE_CHARM */

#endif /* MESH_BLOCK_MPI_HPP */

