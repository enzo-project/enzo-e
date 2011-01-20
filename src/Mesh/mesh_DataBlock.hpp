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
  DataBlock() throw()
    : field_block_(new FieldBlock)
  {  };

  /// Delete DataBlock

  ~DataBlock() throw()
  { delete field_block_; }

  /// Return the Field block

  FieldBlock * field_block ()
  { return field_block_; };

private: // attributes

  FieldBlock    * field_block_;

};

#endif /* MESH_DATA_BLOCK_HPP */

