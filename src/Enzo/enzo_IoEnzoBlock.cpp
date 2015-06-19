// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_IoEnzoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoEnzoBlock class

#include "enzo.hpp"

//----------------------------------------------------------------------

IoEnzoBlock::IoEnzoBlock() throw ()
  : IoBlock(),
    index_enzo_(0)
{
  // save number of meta data elements in IoBlock's

  index_enzo_ = meta_count();

  meta_name_.push_back("enzo_dt");
  meta_name_.push_back("enzo_GridLeftEdge");
  meta_name_.push_back("enzo_GridDimension");
  meta_name_.push_back("enzo_GridStartIndex");
  meta_name_.push_back("enzo_GridEndIndex");
  meta_name_.push_back("enzo_CellWidth");
}

//----------------------------------------------------------------------

void IoEnzoBlock::meta_value
(int index,
 void ** buffer,  std::string * name, scalar_type * type,
 int * nxd, int * nyd, int * nzd) throw()
{

  if (index < index_enzo_) {

    IoBlock::meta_value(index,buffer,name,type,nxd,nyd,nzd);

  } else {

    Io::meta_value(index,buffer,name,type,nxd,nyd,nzd);

    const EnzoBlock * enzo_block = dynamic_cast<const EnzoBlock *>(block_);

    int index_count = index_enzo_;

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
