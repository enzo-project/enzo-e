// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoBlock class

#ifndef IO_IO_BLOCK_HPP
#define IO_IO_BLOCK_HPP

class CommBlock;

class IoBlock : public Io {

  /// @class    IoBlock
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between CommBlock and Output classes

public: // interface

  /// Constructor
  IoBlock() throw();

  /// Destructor
  virtual ~IoBlock() throw()
  {}

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    // NOTE: change this function whenever attributes change
    Io::pup(p);
    WARNING ("IoBlock::pup","skipping block_");
    //    p | *block_;
  }
#endif

  /// Set block
  void set_block (CommBlock * block) throw()
  { block_ = block; };

#include "_io_Io_common.hpp"
  
protected: // functions

protected: // attributes

  CommBlock * block_;

};

#endif /* IO_IO_BLOCK_HPP */

