// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_IoEnzoSimulation.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-05-15
/// @brief    Implementation of the IoEnzoSimulation class

#include "enzo.hpp"

//----------------------------------------------------------------------

IoEnzoSimulation::IoEnzoSimulation(const EnzoSimulation * s) throw ()
  : IoSimulation(s),
    index_enzo_(),
    enzo_turbou_real_state_(s->turbou_real_state_),
    enzo_turbou_int_state_(s->turbou_int_state_),
    num_real_(s->turbou_real_state_.size()),
    num_int_(s->turbou_int_state_.size())
{
  // save number of meta data elements in IoSimulation's

  index_enzo_ = meta_count();

  if (num_real_ > 0) {
    meta_name_.push_back("enzo_turbou_num_real");
    meta_name_.push_back("enzo_turbou_num_int");
    meta_name_.push_back("enzo_turbou_real_state");
    meta_name_.push_back("enzo_turbou_int_state");
  }
}

//----------------------------------------------------------------------

void IoEnzoSimulation::meta_value
(int index,
 void ** buffer,  std::string * name, int * type,
 int * nxd, int * nyd, int * nzd) throw()
{

  if (index < index_enzo_) {

    IoSimulation::meta_value(index,buffer,name,type,nxd,nyd,nzd);

  } else {

    Io::meta_value(index,buffer,name,type,nxd,nyd,nzd);

    int index_count = index_enzo_;

    if (index == index_count++) {

      *buffer = (void *) &num_real_;
      *type   = type_int;
      *nxd    = 1;

    } else if (index == index_count++) {

      *buffer = (void *) &num_int_;
      *type   = type_int;
      *nxd    = 1;

    } else if (index == index_count++) {

      if (enzo_turbou_real_state_.size() < num_real_) {
        enzo_turbou_real_state_.resize(num_real_);
      }
      *buffer = (void *) enzo_turbou_real_state_.data();
      *type   = type_double;
      *nxd     = num_real_;

    } else if (index == index_count++) {
      if (enzo_turbou_int_state_.size() < num_int_) {
        enzo_turbou_int_state_.resize(num_int_);
      }
      *buffer = (void *) enzo_turbou_int_state_.data();
      *type   = type_int;
      *nxd    = num_int_;

    }
  }
}
//----------------------------------------------------------------------

void IoEnzoSimulation::data_value
(int index,
 void ** buffer, std::string * name, int * type,
 int * nxd, int * nyd, int * nzd) throw()
{
}

//======================================================================

int IoEnzoSimulation::data_size () const
{
  int size = 0;

  size += IoSimulation::data_size();

  SIZE_SCALAR_TYPE(size,int,num_real_);
  SIZE_SCALAR_TYPE(size,int,num_int_);
  SIZE_VECTOR_TYPE(size,double,enzo_turbou_real_state_);
  SIZE_VECTOR_TYPE(size,int,enzo_turbou_int_state_);

  return size;
}

//----------------------------------------------------------------------

char * IoEnzoSimulation::save_data (char * buffer) const
{
  char * pc = buffer;

  pc = IoSimulation::save_data(pc);
  
  SAVE_SCALAR_TYPE(pc,int,num_real_);
  SAVE_SCALAR_TYPE(pc,int,num_int_);
  SAVE_VECTOR_TYPE(pc,double,enzo_turbou_real_state_);
  SAVE_VECTOR_TYPE(pc,int,enzo_turbou_int_state_);

  ASSERT2 ("IoEnzoSimulation::save_data()",
  	   "Expecting buffer size %d actual size %d",
  	   IoEnzoSimulation::data_size(),(pc-buffer),
  	   (IoEnzoSimulation::data_size() == (pc-buffer)));
  
  // return first byte after filled buffer
  return pc;
}

//----------------------------------------------------------------------

char * IoEnzoSimulation::load_data (char * buffer)
{
  char * pc = buffer;

  pc = IoSimulation::load_data(pc);
  
  LOAD_SCALAR_TYPE(pc,int,num_real_);
  LOAD_SCALAR_TYPE(pc,int,num_int_);
  LOAD_VECTOR_TYPE(pc,double,enzo_turbou_real_state_);
  LOAD_VECTOR_TYPE(pc,int,enzo_turbou_int_state_);

  return pc;
}

//----------------------------------------------------------------------

void IoEnzoSimulation::save_to (void * v)
{
  IoSimulation::save_to(v);
  
  EnzoSimulation * enzo_simulation = (EnzoSimulation *)v;

  enzo_simulation->turbou_real_state_ = enzo_turbou_real_state_;
  enzo_simulation->turbou_int_state_ = enzo_turbou_int_state_;

}

//----------------------------------------------------------------------
