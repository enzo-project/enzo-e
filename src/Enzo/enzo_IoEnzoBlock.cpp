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
  meta_name_.push_back("enzo_redshift");
  
}

void IoEnzoBlock::set_block (Block * block) throw()
{

  IoBlock::set_block (block);

  EnzoBlock * enzo_block = static_cast<EnzoBlock*>(block);
  
  enzo_dt_ = enzo_block->dt;
  for (int i=0; i<3; i++) {
    enzo_GridLeftEdge_[i] = enzo_block->GridLeftEdge[i];
    enzo_GridDimension_[i] = enzo_block->GridDimension[i];
    enzo_GridStartIndex_[i] = enzo_block->GridStartIndex[i];
    enzo_GridEndIndex_[i] = enzo_block->GridEndIndex[i];
    enzo_CellWidth_[i] = enzo_block->CellWidth[i];
  }
  enzo_redshift_ = enzo_block->redshift;
    
}

//----------------------------------------------------------------------

void IoEnzoBlock::meta_value
(int index,
 void ** buffer,  std::string * name, int * type,
 int * nxd, int * nyd, int * nzd) throw()
{

  if (index < index_enzo_) {

    IoBlock::meta_value(index,buffer,name,type,nxd,nyd,nzd);

  } else {

    Io::meta_value(index,buffer,name,type,nxd,nyd,nzd);

    int index_count = index_enzo_;

    if (index == index_count++) {

      *buffer = (void *) & enzo_dt_;
      *type   = type_enzo_float;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_GridLeftEdge_;
      *type   = type_enzo_float;
      *nxd     = 3;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_GridDimension_;
      *type   = type_int;
      *nxd     = 3;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_GridStartIndex_;
      *type   = type_int;
      *nxd     = 3;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_GridEndIndex_;
      *type   = type_int;
      *nxd     = 3;

    } else if (index == index_count++) {

      *buffer = (void *) enzo_CellWidth_;
      *type   = type_enzo_float;
      *nxd     = 3;
      
    } else if (index == index_count++) {
      *buffer = (void *) & enzo_redshift_;
      *type   = type_enzo_float;
    }
  }
}
//----------------------------------------------------------------------

void IoEnzoBlock::data_value
(int index,
 void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd) throw()
{
}

//======================================================================

int IoEnzoBlock::data_size () const
{
  int size = 0;

  size += IoBlock::data_size();

  SIZE_SCALAR_TYPE(size,int, index_enzo_);
  SIZE_SCALAR_TYPE(size,enzo_float,enzo_dt_);
  SIZE_ARRAY_TYPE(size,enzo_float,enzo_GridLeftEdge_,3);
  SIZE_ARRAY_TYPE(size,int, enzo_GridDimension_,3);
  SIZE_ARRAY_TYPE(size,int, enzo_GridStartIndex_,3);
  SIZE_ARRAY_TYPE(size,int, enzo_GridEndIndex_,3);
  SIZE_ARRAY_TYPE(size,enzo_float,enzo_CellWidth_,3);
  SIZE_SCALAR_TYPE(size,enzo_float,enzo_redshift_);

  return size;
}

//----------------------------------------------------------------------

char * IoEnzoBlock::save_data (char * buffer) const
{
  char * pc = buffer;

  pc = IoBlock::save_data(pc);
  
  SAVE_SCALAR_TYPE(pc,int, index_enzo_);
  SAVE_SCALAR_TYPE(pc,enzo_float,enzo_dt_);
  SAVE_ARRAY_TYPE(pc,enzo_float,enzo_GridLeftEdge_,3);
  SAVE_ARRAY_TYPE(pc,int, enzo_GridDimension_,3);
  SAVE_ARRAY_TYPE(pc,int, enzo_GridStartIndex_,3);
  SAVE_ARRAY_TYPE(pc,int, enzo_GridEndIndex_,3);
  SAVE_ARRAY_TYPE(pc,enzo_float,enzo_CellWidth_,3);
  SAVE_SCALAR_TYPE(pc,enzo_float,enzo_redshift_);

  ASSERT2 ("IoEnzoBlock::save_data()",
  	   "Expecting buffer size %d actual size %d",
  	   IoEnzoBlock::data_size(),(pc-buffer),
  	   (IoEnzoBlock::data_size() == (pc-buffer)));
  
  // return first byte after filled buffer
  return pc;
}

//----------------------------------------------------------------------

char * IoEnzoBlock::load_data (char * buffer)
{
  char * pc = buffer;

  pc = IoBlock::load_data(pc);
  
  LOAD_SCALAR_TYPE(pc,int, index_enzo_);
  LOAD_SCALAR_TYPE(pc,enzo_float,enzo_dt_);
  LOAD_ARRAY_TYPE(pc,enzo_float,enzo_GridLeftEdge_,3);
  LOAD_ARRAY_TYPE(pc,int, enzo_GridDimension_,3);
  LOAD_ARRAY_TYPE(pc,int, enzo_GridStartIndex_,3);
  LOAD_ARRAY_TYPE(pc,int, enzo_GridEndIndex_,3);
  LOAD_ARRAY_TYPE(pc,enzo_float,enzo_CellWidth_,3);
  LOAD_SCALAR_TYPE(pc,enzo_float,enzo_redshift_);

  return pc;
}

//----------------------------------------------------------------------
