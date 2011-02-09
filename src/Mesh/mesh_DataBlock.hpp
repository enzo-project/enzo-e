// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_DataBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Mesh] Declaration of the DataBlock class

#ifndef MESH_DATA_BLOCK_HPP
#define MESH_DATA_BLOCK_HPP

class DataBlock {

  /// @class    DataBlock
  /// @ingroup  Data
  /// @brief    [\ref Mesh] Container class for all data blocks (currently just fields)

public: // interface

  /// Initialize the DataBlock object
  DataBlock(int num_field_blocks = 1) throw()
    : field_block_()
  { field_block_.resize(num_field_blocks);
    for (size_t i=0; i<field_block_.size(); i++) {
      field_block_[i] = new FieldBlock;
    }
  }

  /// Delete DataBlock

  ~DataBlock() throw()
  { 
    for (size_t i=0; i<field_block_.size(); i++) {
      delete field_block_[i];
    }
  }

  /// Return the Field block

  FieldBlock * field_block (int i=0)
  { return field_block_.at(i); };

private: // attributes

  /// Array of field blocks

  std::vector<FieldBlock *> field_block_;

};

#endif /* MESH_DATA_BLOCK_HPP */

