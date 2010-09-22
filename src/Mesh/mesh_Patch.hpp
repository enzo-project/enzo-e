// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_PATCH_HPP
#define MESH_PATCH_HPP

/// @file     mesh_Patch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the interface for the Patch class

class Patch {

  /// @class    Patch
  /// @ingroup  Mesh
  /// @brief    Brief description of class Patch.

public: // interface

  /// Constructor
  Patch(DataDescr * data_descr,
	int patch_size[3], 
	int block_size,
	double lower[3],
	double upper[3]) ;

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Patch() throw();

  /// Copy constructor
  Patch(const Patch & patch) throw();

  /// Assignment operator
  Patch & operator= (const Patch & patch) throw();

  //--------------------------------------------------

  /// Initial patch advance, ending with receive_()
  void p_evolve();

  /// Return the number of local data blocks
  int num_data_blocks()  throw()
  { return block_count_[0]*block_count_[1]*block_count_[2]; };

  /// Return the ith local data block
  DataBlock * data_block(int ix, int iy, int iz) throw()
  { return data_block_[ix + block_count_[0] * (iy + block_count_[1] * iz)]; };

private: // attributes

  /// Array of local data blocks associated with the patch
  DataBlock ** data_block_;

  /// Size of the patch
  int patch_size_[3];

  /// Number of blocks
  int block_count_[3];

  /// Size of each block
  int block_size_;

  /// Lower extent of the patch
  double lower_[3];

  /// Upper extent of the patch
  double upper_[3];

};

#endif /* MESH_PATCH_HPP */

