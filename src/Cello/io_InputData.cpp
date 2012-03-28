// See LICENSE_CELLO file for license and copyright information

/// @file     io_InputData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-23
/// @brief    Implementation of the InputData class

#include "cello.hpp"
#include "io.hpp"

//----------------------------------------------------------------------

InputData::InputData(const Factory * factory) throw ()
  : Input(factory)
{
}

//======================================================================

void InputData::open () throw()
{
  std::string file_name = expand_file_name_(&file_name_,&file_args_);

  Monitor::instance()->print ("Input","reading data file %s", 
			      file_name.c_str());

  close();

  file_ = new FileHdf5 (".",file_name);

  file_->file_open();
}

//----------------------------------------------------------------------

InputData::~InputData() throw()
{
  close();
}

//----------------------------------------------------------------------

bool InputData::is_open () throw()
{
  return (file_ != 0);
}

//----------------------------------------------------------------------

void InputData::close () throw()
{
  if (file_) file_->file_close();

  delete file_;  file_ = 0;
}

//----------------------------------------------------------------------

void InputData::finalize () throw ()
{
  Input::finalize();
}

// //----------------------------------------------------------------------

// void InputData::read_hierarchy 
// (
//  Hierarchy  * hierarchy,
//  const FieldDescr * field_descr
//  ) throw()
// {

//   IoHierarchy io_hierarchy(hierarchy);

//   // Read hierarchy meta-data

  
//   // file_->file_read_meta("value", "name",scalar_type_char,6);

//   Input::read_meta (&io_hierarchy);

//   // Call read_patch() on contained patches
//   Input::read_hierarchy (hierarchy, field_descr);

// }

//----------------------------------------------------------------------

void InputData::read_patch 
(
 Patch * patch,
 const FieldDescr * field_descr,
 int ixp0, int iyp0, int izp0
 ) throw()
{
  // Create file group for patch

  char buffer[40];
  int ib = patch->index();
  sprintf (buffer,"patch_%d",ib);
  file_->group_chdir(buffer);
  file_->group_create();

  // Loop over metadata items in Hierarchy

  IoPatch io_patch(patch);

  // Read patch meta-data
  Input::read_meta_group (&io_patch);

  // Call read_block() on contained blocks
  Input::read_patch(patch,field_descr,ixp0,iyp0,izp0);

#ifndef CONFIG_USE_CHARM
  file_->group_close();
  file_->group_chdir("..");
#endif

}

#ifdef CONFIG_USE_CHARM
//----------------------------------------------------------------------

void InputData::end_read_patch() throw()
{
  file_->group_close();
  file_->group_chdir("..");
}
//----------------------------------------------------------------------

#endif

void InputData::read_block 
( 
 Block * block,
  const FieldDescr * field_descr,
  int ixp0, int iyp0, int izp0) throw()
{

  // Create file group for block

  char buffer[40];
  int ib = block->index();
  sprintf (buffer,"block_%d",ib);
  file_->group_chdir(buffer);
  file_->group_create();

  // Read block meta data

  io_block()->set_block(block);

  Input::read_meta_group (io_block());

  // Call read_block() on base Input object

  Input::read_block(block,field_descr,ixp0,iyp0,izp0);

  file_->group_close();
  file_->group_chdir("..");

}

//----------------------------------------------------------------------

void InputData::read_field
( 
 FieldBlock * field_block,
 const FieldDescr * field_descr,
 int field_index) throw()
{
  io_field_block()->set_field_descr(field_descr);
  io_field_block()->set_field_block(field_block);
  io_field_block()->set_field_index(field_index);

  for (size_t i=0; i<io_field_block()->data_count(); i++) {

    void * buffer;
    std::string name;
    scalar_type type;
    int nxd,nyd,nzd;  // Array dimension
    int nx,ny,nz;     // Array size

    // Get ith FieldBlock data
    io_field_block()->data_value(i, &buffer, &name, &type, 
				 &nxd,&nyd,&nzd,
				 &nx, &ny, &nz);

    // Read ith FieldBlock data

    file_->data_create(name.c_str(),type,nxd,nyd,nzd,nx,ny,nz);
    file_->data_read(buffer);
    file_->data_close();
  }

}

//======================================================================
