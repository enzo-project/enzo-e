// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_IoEnzoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoEnzoBlock class

#include "enzo.hpp"

//----------------------------------------------------------------------

IoEnzoBlock::IoEnzoBlock() throw ()
  : IoBlock()
{

  meta_name_.push_back("enzo_dt");
  meta_count_ ++;
  meta_name_.push_back("enzo_GridLeftEdge");
  meta_count_ ++;
  meta_name_.push_back("enzo_GridDimension");
  meta_count_ ++;
  meta_name_.push_back("enzo_GridStartIndex");
  meta_count_ ++;
  meta_name_.push_back("enzo_GridEndIndex");
  meta_count_ ++;
  meta_name_.push_back("enzo_CellWidth");
  meta_count_ ++;
  
  ASSERT2("IoEnzoBlock::IoEnzoBlock()",
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
    int index_count = 0;

    if (index == index_count++) {

      *buffer = (void *) & enzo_block->dt;
      *type   = scalar_type_enzo_float;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_block->GridLeftEdge;
      *type   = scalar_type_enzo_float;
      *nxd     = 3;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_block->GridDimension;
      *type   = scalar_type_int;
      *nxd     = 3;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_block->GridStartIndex;
      *type   = scalar_type_int;
      *nxd     = 3;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_block->GridEndIndex;
      *type   = scalar_type_int;
      *nxd     = 3;

    } else if (index == index_count++) {

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
