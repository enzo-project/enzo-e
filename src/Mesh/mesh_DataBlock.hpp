// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_DataBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @todo     Move FieldBlock members//attributes to DataBlock when possible
/// @brief    [\ref Mesh] Declaration of the DataBlock class

#ifndef MESH_DATA_BLOCK_HPP
#define MESH_DATA_BLOCK_HPP

class DataBlock {

  /// @class    DataBlock
  /// @ingroup  Data
  /// @brief    [\ref Mesh] Container class for all data blocks (currently just fields)

public: // interface

  /// Initialize the DataBlock object
  DataBlock(DataDescr * data_descr,
	    int num_field_blocks = 1) throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~DataBlock() throw();

  /// Copy constructor
  DataBlock(const DataBlock & data_block) throw();

  /// Assignment operator
  DataBlock & operator= (const DataBlock & data_block) throw();

  //----------------------------------------------------------------------

  /// Return the Field block
  const FieldBlock * field_block (int i=0) const throw()
  { return field_block_.at(i); };

  /// Return the Field block
  FieldBlock * field_block (int i=0) throw()
  { return field_block_.at(i); };

private: // functions

  /// Allocate and copy in attributes from data_block
  void create_(const DataBlock & data_block) throw();

private: // attributes

  /// Pointer to the parent DataDescr
  DataDescr * data_descr_;

  /// Array of field blocks
  std::vector<FieldBlock *> field_block_;

};

#endif /* MESH_DATA_BLOCK_HPP */

