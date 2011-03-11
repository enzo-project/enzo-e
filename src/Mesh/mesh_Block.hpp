// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @todo     Move FieldBlock members//attributes to Block when possible
/// @todo     Add Patch pointer, and exchang computable data for function calls
/// @brief    [\ref Mesh] Declaration of the Block class

#ifndef MESH_BLOCK_HPP
#define MESH_BLOCK_HPP

class FieldDescr;
class FieldBlock;
class Patch;

class Block {

  /// @class    Block
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Basic serial block of mesh data

public: // interface

  /// Initialize the Block object
  Block(Patch * patch,
	FieldDescr * field_descr,
	int nx, int ny=1, int nz=1,
	int num_field_blocks = 1) throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  virtual ~Block() throw();

  /// Copy constructor
  Block(const Block & block) throw();

  /// Assignment operator
  Block & operator= (const Block & block) throw();

  //----------------------------------------------------------------------

  /// Return the ith Field block
  const FieldBlock * field_block (int i=0) const throw();

  /// Return the ith Field block
  FieldBlock * field_block (int i=0) throw();

  /// Return lower values of the block (excluding ghosts)
  void extent(double * lower_x = 0, double * upper_x = 0, 
	      double * lower_y = 0, double * upper_y = 0,
	      double * lower_z = 0, double * upper_z = 0) const throw ();

  /// Set the box extent
  void set_extent(double lower_x = 0.0, double upper_x = 1.0, 
		  double lower_y = 0.0, double upper_y = 1.0,
		  double lower_z = 0.0, double upper_z = 1.0) throw();

protected: // functions

  /// Allocate and copy in attributes from give Block
  void copy_(const Block & block) throw();

protected: // attributes

  /// Parent Patch
  Patch * patch_;
 
  /// Array of field blocks
  std::vector<FieldBlock *> field_block_;

  /// Extent of the box associated with the block
  /// WARNING: should not be used for deep AMR due to precision /
  /// range issues
  double lower_[3];
  double upper_[3];

};

#endif /* MESH_BLOCK_HPP */

