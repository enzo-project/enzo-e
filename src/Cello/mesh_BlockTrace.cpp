// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_BlockTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-26
/// @brief    Implementation of the BlockTrace class

#include "mesh.hpp"

//----------------------------------------------------------------------

int BlockTrace::data_size () const
{
  int size = 0;
  SIZE_SCALAR_TYPE(size,int,num_children_);
  SIZE_SCALAR_TYPE(size,Index,index_home_);
  SIZE_SCALAR_TYPE(size,Index,index_root_);
  SIZE_ARRAY_TYPE(size,int,index_min_,3);
  SIZE_ARRAY_TYPE(size,int,index_max_,3);
  SIZE_ARRAY_TYPE(size,int,index_curr_,3);
  SIZE_VECTOR_TYPE(size,int,child_stack_);
  SIZE_VECTOR_TYPE(size,Index,index_stack_);
  return size;
}

//----------------------------------------------------------------------

char * BlockTrace::save_data (char * buffer) const
{
  char * pc = buffer;
  SAVE_SCALAR_TYPE(pc,int,num_children_);
  SAVE_SCALAR_TYPE(pc,Index,index_home_);
  SAVE_SCALAR_TYPE(pc,Index,index_root_);
  SAVE_ARRAY_TYPE(pc,int,index_min_,3);
  SAVE_ARRAY_TYPE(pc,int,index_max_,3);
  SAVE_ARRAY_TYPE(pc,int,index_curr_,3);
  SAVE_VECTOR_TYPE(pc,int,child_stack_);
  SAVE_VECTOR_TYPE(pc,Index,index_stack_);
  ASSERT2 ("DataMsg::save_data()",
           "Expecting buffer size %d actual size %d",
           data_size(),(pc-buffer),
           (data_size() == (pc-buffer)));

  return pc;
}

//----------------------------------------------------------------------

char * BlockTrace::load_data (char * buffer)
{
  char * pc = buffer;
  LOAD_SCALAR_TYPE(pc,int,num_children_);
  LOAD_SCALAR_TYPE(pc,Index,index_home_);
  LOAD_SCALAR_TYPE(pc,Index,index_root_);
  LOAD_ARRAY_TYPE(pc,int,index_min_,3);
  LOAD_ARRAY_TYPE(pc,int,index_max_,3);
  LOAD_ARRAY_TYPE(pc,int,index_curr_,3);
  LOAD_VECTOR_TYPE(pc,int,child_stack_);
  LOAD_VECTOR_TYPE(pc,Index,index_stack_);
  ASSERT2 ("DataMsg::load_data()",
           "Expecting buffer size %d actual size %d",
           data_size(),(pc-buffer),
           (data_size() == (pc-buffer)));
  return pc;
}
