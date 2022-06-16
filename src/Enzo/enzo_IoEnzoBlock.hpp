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

  /// CHARM++ PUP::able declaration
  PUPable_decl(IoEnzoBlock);

  /// CHARM++ migration constructor
  IoEnzoBlock(CkMigrateMessage *m) : IoBlock(m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    IoBlock::pup(p);
    TRACEPUP;

    p | index_enzo_;

  }

  /// Return the ith metadata item associated with the EnzoBlock object
  void meta_value 
  (int index, 
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Return the ith data item associated with the EnzoBlock object
  void data_value 
  (int index, 
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  virtual void set_block (Block * block) throw();

  /// PACKING / UNPACKING
  
  /// Return the number of bytes required to serialize the data object
  virtual int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  virtual char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  virtual char * load_data (char * buffer);

  virtual void print(const char * msg)
  {
    CkPrintf ("DEBUG_IO_ENZO_BLOCK %d %s\n",CkMyPe(),msg);
    CkPrintf ("DEBUG_IO_ENZO_BLOCK index_enzo_          %d\n", index_enzo_);
    CkPrintf ("DEBUG_IO_ENZO_BLOCK enzo_dt_             %g\n", enzo_dt_);
    CkPrintf ("DEBUG_IO_ENZO_BLOCK enzo_GridLeftEdge_   %g %g %g\n",
              enzo_GridLeftEdge_[0], enzo_GridLeftEdge_[1], enzo_GridLeftEdge_[2]);
    CkPrintf ("DEBUG_IO_ENZO_BLOCK enzo_GridDimension_  %d %d %d\n",
              enzo_GridDimension_[0], enzo_GridDimension_[1], enzo_GridDimension_[2]);
    CkPrintf ("DEBUG_IO_ENZO_BLOCK enzo_GridStartIndex_ %d %d %d\n",
              enzo_GridStartIndex_[0], enzo_GridStartIndex_[1], enzo_GridStartIndex_[2]);
    CkPrintf ("DEBUG_IO_ENZO_BLOCK enzo_GridEndIndex_   %d %d %d\n",
              enzo_GridEndIndex_[0], enzo_GridEndIndex_[1], enzo_GridEndIndex_[2]);
    CkPrintf ("DEBUG_IO_ENZO_BLOCK enzo_CellWidth_      %g %g %g\n",
              enzo_CellWidth_[0], enzo_CellWidth_[1], enzo_CellWidth_[2]);
    CkPrintf ("DEBUG_IO_ENZO_BLOCK eonz_redshift_       %g\n",
              enzo_redshift_);

    IoBlock::print(msg);
  }
private:

  int index_enzo_;

  enzo_float enzo_dt_;
  enzo_float enzo_GridLeftEdge_[3];
  int enzo_GridDimension_[3];
  int enzo_GridStartIndex_[3];
  int enzo_GridEndIndex_[3];
  enzo_float enzo_CellWidth_[3];
  enzo_float enzo_redshift_;
};

#endif /* IO_IO_ENZO_BLOCK_HPP */

