// $Id: data_DataBlock.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    Declaration of the DataBlock class

#ifndef DATA_DATA_BLOCK_HPP
#define DATA_DATA_BLOCK_HPP

#include <vector>

#include "error.hpp"
#include "field.hpp"

class DataBlock {

  /// @class    DataBlock
  /// @ingroup  Data
  /// @brief    Container class for all data blocks (currently just fields)

public: // interface

  /// Initialize the DataBlock object
  DataBlock(FieldBlock * field_block) throw()
    : field_block_(field_block)
  {
    ASSERT("DataBlock()","field_block must be non-null",field_block != NULL);
  };

  //----------------------------------------------------------------------
  // Field functions
  //----------------------------------------------------------------------

  /// Return the Field block

  FieldBlock * field_block ()
  { return field_block_; };

private: // attributes

  FieldBlock    * field_block_;

};

#endif /* DATA_DATA_BLOCK_HPP */

