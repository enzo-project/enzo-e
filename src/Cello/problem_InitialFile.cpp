// See LICENSE_CELLO file for license and copyright information

/// @file     method_InitialFile.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-02-16
/// @brief    Implementation of the InitialFile class

#include "cello.hpp"

#include "problem.hpp"

//----------------------------------------------------------------------

InitialFile::InitialFile
(Parameters * parameters,
 GroupProcess * group_process,
 int cycle, double time) throw ()
  : Initial (cycle,time),
    parameters_(parameters),
    group_process_(group_process),
    input_(0)
{
  TRACE("InitialFile::InitialFile");
}

//----------------------------------------------------------------------

InitialFile::~InitialFile() throw()
{
  delete input_; input_ = 0;
}

//----------------------------------------------------------------------

void InitialFile::enforce
(
 const Hierarchy  * hierarchy,
 const FieldDescr * field_descr,
 Block            * block
 ) throw()
{
  if (! input_) input_ = new InputData(hierarchy->factory());

  if (! input_->is_open() ) {

    input_->is_scheduled (cycle_,time_);

    std::string              file_name = "";
    std::vector<std::string> file_args;

    get_filename_(&file_name,&file_args);

    input_->open();

  }

  INCOMPLETE("InitialFile::enforce");

  Patch * patch = hierarchy->factory()->create_patch
    (group_process_,
     0,0,0,
     0,0,0,
     0,0,0,
     0,0,0,
     0,0,0);

  input_->read_patch(patch,field_descr,0,0,0);

  // DEBUG

  int nx=0,ny=0,nz=0;
  patch->size(&nx,&ny,&nz);
  TRACE3 ("Patch size %d %d %d",nx,ny,nz);

  int nbx=0,nby=0,nbz=0;
  patch->layout()->block_count(&nbx,&nby,&nbz);
  int p0=0,np=0;
  patch->layout()->process_range(&p0,&np);
  TRACE5 ("Patch layout %d:%d (%d %d %d)",p0,p0+np-1,
	  nbx,nby,nbz);


  nx=0,ny=0,nz=0;
  patch->offset(&nx,&ny,&nz);
  TRACE3 ("Patch offset %d %d %d",nx,ny,nz);

  nx=0,ny=0,nz=0;
  patch->blocking(&nx,&ny,&nz);
  TRACE3 ("Patch blocking %d %d %d",nx,ny,nz);

  double x=0,y=0,z=0;
  patch->lower(&x,&y,&z);
  TRACE3 ("Patch lower %f %f %f",x,y,z);

  x=0,y=0,z=0;
  patch->upper(&x,&y,&z);
  TRACE3 ("Patch upper %f %f %f",x,y,z);

  // (X) Layout * layout_;
  // (*) int size_[3];
  // ( ) int offset_[3];
  // ( ) int blocking_[3];
  // ( ) double lower_[3];
  // ( ) double upper_[3];
}

//----------------------------------------------------------------------

void InitialFile::get_filename_
(
 std::string * file_name,
 std::vector<std::string> * file_args
 ) throw()
{
  // parameter: Initial : name

  parameters_->group_set(0,"Initial");

  if (parameters_->type("name") == parameter_string) {

    *file_name = parameters_->value_string("name","");
     
  } else if (parameters_->type("name") == parameter_list) {

    int list_length = parameters_->list_length("name");

    *file_name = parameters_->list_value_string(0,"name","");

    for (int index = 1; index<list_length; index++) {
      file_args->push_back(parameters_->list_value_string(index,"name",""));
    }

  } else {

    ERROR1("InitialFile::enforce",
	   "Bad type %d for 'Initial : name' parameter",
	   parameters_->type("name"));

  }

  TRACE("Setting file_name,file_args");

  input_->set_filename (*file_name, *file_args);

}
