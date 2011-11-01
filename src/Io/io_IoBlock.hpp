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

  /// Set block
  void set_block (const Block * block) throw()
  { block_ = block; };

  /// Return the ith metadata item associated with the Block object
  void meta_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  /// Return the ith data item associated with the Block object
  void data_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  
protected: // functions

  const Block * block_;

private: // attributes


};

#endif /* IO_IO_BLOCK_HPP */

