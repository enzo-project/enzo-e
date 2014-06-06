// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_IoEnzoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Enzo] Declaration of the IoEnzoBlock class

#ifndef IO_IO_ENZO_BLOCK_HPP
#define IO_IO_ENZO_BLOCK_HPP

class EnzoBlock;

class IoEnzoBlock : public IoBlock {

  /// @class    IoEnzoBlock
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Class for interfacing between EnzoBlock and Output objects

public: // interface

  /// Constructor
  IoEnzoBlock() throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change

    TRACEPUP;

    IoBlock::pup(p);

    p | meta_count_enzo_;

  }

  /// Return the ith metadata item associated with the EnzoBlock object
  void meta_value 
  (int index, 
   void ** buffer, std::string * name, scalar_type * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Return the ith data item associated with the EnzoBlock object
  void data_value 
  (int index, 
   void ** buffer, std::string * name, scalar_type * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

private:
  int meta_count_enzo_;

};

#endif /* IO_IO_ENZO_BLOCK_HPP */

