// See LICENSE_CELLO file for license and copyright information

/// @file     Io_IoEnzoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoEnzoBlock class

#include "enzo.hpp"

//----------------------------------------------------------------------

IoEnzoBlock::IoEnzoBlock() throw ()
  : IoBlock()
{
  meta_count_enzo_ = 9;

  meta_count_ += meta_count_enzo_;

  meta_name_.push_back("enzo_Time");
  meta_name_.push_back("enzo_CycleNumber");
  meta_name_.push_back("enzo_OldTime");
  meta_name_.push_back("enzo_dt");
  meta_name_.push_back("enzo_GridLeftEdge");
  meta_name_.push_back("enzo_GridDimension");
  meta_name_.push_back("enzo_GridStartIndex");
  meta_name_.push_back("enzo_GridEndIndex");
  meta_name_.push_back("enzo_CellWidth");
  
  ASSERT2("IoPatch::IoPatch()",
	 "meta_name.size() [%d] !=  meta_count_ [%d]",
	  meta_name_.size(),meta_count(),
	  meta_name_.size()==meta_count());

}

//----------------------------------------------------------------------

void IoEnzoBlock::meta_value
(int index,
 void ** buffer,  std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{

  int index_block = meta_count_ - meta_count_enzo_;

  if (index < index_block) {

    // First return Block attributes
    IoBlock::meta_value(index,buffer,name,type,nxd,nyd,nzd);

  } else {

    // Then return EnzoBlock attributes

    Io::meta_value(index,buffer,name,type,nxd,nyd,nzd);

    const EnzoBlock * enzo_block = dynamic_cast<const EnzoBlock *>(block_);

    index -= index_block;

    if (index == 0) {

      *buffer = (void *) & enzo_block->Time_;
      *type   = scalar_type_enzo_float;

    } else if (index == 1) {

      *buffer = (void *) & enzo_block->CycleNumber;
      *type   = scalar_type_int;

    } else if (index == 2) {

      *buffer = (void *) & enzo_block->OldTime;
      *type   = scalar_type_enzo_float;

    } else if (index == 3) {

      *buffer = (void *) & enzo_block->dt;
      *type   = scalar_type_enzo_float;

    } else if (index == 4) {

      *buffer = (void *) enzo_block->GridLeftEdge;
      *type   = scalar_type_enzo_float;
      *nxd     = 3;

    } else if (index == 5) {

      *buffer = (void *) enzo_block->GridDimension;
      *type   = scalar_type_int;
      *nxd     = 3;

    } else if (index == 6) {

      *buffer = (void *) enzo_block->GridStartIndex;
      *type   = scalar_type_int;
      *nxd     = 3;

    } else if (index == 7) {

      *buffer = (void *) enzo_block->GridEndIndex;
      *type   = scalar_type_int;
      *nxd     = 3;

    } else if (index == 8) {

      *buffer = (void *) enzo_block->CellWidth;
      *type   = scalar_type_enzo_float;
      *nxd     = 3;

    }
  }
}
//----------------------------------------------------------------------

void IoEnzoBlock::data_value
(int index,
 void ** buffer, std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{
}

//----------------------------------------------------------------------
