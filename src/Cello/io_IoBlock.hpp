// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoBlock class

#ifndef IO_IO_BLOCK_HPP
#define IO_IO_BLOCK_HPP

class Block;

class IoBlock : public Io {

  /// @class    IoBlock
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between Block and Output classes

public: // interface

  /// Constructor
  IoBlock() throw();

  /// Destructor
  virtual ~IoBlock() throw()
  {}

  /// Set block
  void set_block (const Block * block) throw()
  { block_ = block; };

#include "_io_Io_common.hpp"
  
protected: // functions

  const Block * block_;

private: // attributes


};

#endif /* IO_IO_BLOCK_HPP */

