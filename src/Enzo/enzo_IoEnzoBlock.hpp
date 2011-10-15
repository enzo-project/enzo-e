// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoEnzoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoEnzoBlock class

#ifndef IO_IO_ENZO_BLOCK_HPP
#define IO_IO_ENZO_BLOCK_HPP

class EnzoBlock;

class IoEnzoBlock : public IoBlock {

  /// @class    IoEnzoBlock
  /// @ingroup  Io
  /// @brief    [\ref Io] 

public: // interface

  /// Constructor
  IoEnzoBlock() throw();

  /// Return the ith metadata item associated with the EnzoBlock object
  void meta_value 
  (int index, 
   void ** buffer, const char ** name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  /// Return the ith data item associated with the EnzoBlock object
  void data_value 
  (int index, 
   void ** buffer, const char ** name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

private:
  int meta_count_enzo_;

};

#endif /* IO_IO_ENZO_BLOCK_HPP */

