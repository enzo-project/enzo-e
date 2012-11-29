// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-28
/// @brief    [\ref Comm] Declaration of the CommBlock class

#ifndef COMM_COMM_BLOCK_HPP
#define COMM_COMM_BLOCK_HPP

class CommBlock {

  /// @class    CommBlock
  /// @ingroup  Comm
  /// @brief    [\ref Comm] 

public: // interface

  /// Constructor
  CommBlock() throw();

  /// Destructor
  ~CommBlock() throw();

  /// Copy constructor
  CommBlock(const CommBlock & CommBlock) throw();

  /// Assignment operator
  CommBlock & operator= (const CommBlock & CommBlock) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
  }
#endif
  
private: // functions


private: // attributes

  /// Block that this CommBlock is associated with
  Block * block_;
  // NOTE: change pup() function whenever attributes change

};

#endif /* COMM_COMM_BLOCK_HPP */

