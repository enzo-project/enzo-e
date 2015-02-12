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

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    // NOTE: change this function whenever attributes change
    Io::pup(p);
    WARNING ("IoBlock::pup","skipping block_");
    //    p | *block_;
  }

  /// Set block
  void set_block (Block * block) throw()
  { block_ = block; };

#include "_io_Io_common.hpp"
  
protected: // functions

protected: // attributes

  Block * block_;

};

#endif /* IO_IO_BLOCK_HPP */

