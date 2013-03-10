// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-03-10
/// @brief    [\ref Mesh] Declaration of the Block class
///

#ifndef MESH_BLOCK_HPP
#define MESH_BLOCK_HPP

class Block {

  /// @class    Block
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Block() throw();

  /// Destructor
  ~Block() throw();

  /// Copy constructor
  Block(const Block & block) throw();

  /// Assignment operator
  Block & operator= (const Block & block) throw();

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

  // NOTE: change pup() function whenever attributes change

};

#endif /* MESH_BLOCK_HPP */

