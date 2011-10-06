// See LICENSE_CELLO file for license and copyright information

/// @file     Io_IoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    Implementation of the IoBlock class

#include "io.hpp"

//----------------------------------------------------------------------

IoBlock::IoBlock(const Block * block) throw ()
  : Io(7,0),
    block_(block)
{
  meta_name_.push_back("size");
  meta_name_.push_back("lower");
  meta_name_.push_back("upper");
  meta_name_.push_back("cycle");
  meta_name_.push_back("time");
  meta_name_.push_back("dt");
#ifdef CONFIG_USE_CHARM
  meta_name_.push_back("index");
#else
  meta_name_.push_back("num_field_blocks");
#endif
}


//----------------------------------------------------------------------

void IoBlock::meta_value
(int index,
 void ** buffer, const char ** name, enum scalar_type * type,
 int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  Io::meta_value(index,buffer,name,type,n0,n1,n2,n3,n4);

  if (index == 0) {

    *buffer = (void *) block_->size_;
    *type   = scalar_type_int;
    *n0     = 3;
    
  } else if (index == 1) {

    *buffer = (void *) block_->lower_;
    *type   = scalar_type_double;
    *n0     = 3;

  } else if (index == 2) {

    *buffer = (void *) block_->upper_;
    *type   = scalar_type_double;
    *n0     = 3;

  } else if (index == 3) {

    *buffer = (void *) & block_->cycle_;
    *type   = scalar_type_int;

  } else if (index == 4) {

    *buffer = (void *) & block_->time_;
    *type   = scalar_type_double;

  } else if (index == 5) {

    *buffer = (void *) & block_->dt_;
    *type   = scalar_type_double;

  } else if (index == 6) {

#ifdef CONFIG_USE_CHARM

    *buffer = (void *) & block_->num_field_blocks_;
    *type   = scalar_type_int;

#else 

    *buffer = (void *) & block_->index_;
    *type   = scalar_type_int;

#endif

  }
}

//----------------------------------------------------------------------

void IoBlock::data_value
(int index,
 void ** buffer, const char ** name, enum scalar_type * type,
 int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
}

//----------------------------------------------------------------------
