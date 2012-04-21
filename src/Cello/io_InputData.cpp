// See LICENSE_CELLO file for license and copyright information

/// @file     io_InputData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-23
/// @brief    Implementation of the InputData class

#include "cello.hpp"
#include "io.hpp"

//----------------------------------------------------------------------

InputData::InputData(const Factory * factory) throw ()
  : Input(factory),
    factory_(factory)
{
}

//----------------------------------------------------------------------

InputData::~InputData() throw()
{
  close();
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

void InputData::close () throw()
{
  if (file_) file_->file_close();

  delete file_;  file_ = 0;
}

//----------------------------------------------------------------------

bool InputData::is_open () throw()
{
  return (file_ != 0);
}

//----------------------------------------------------------------------

void InputData::finalize () throw ()
{
  Input::finalize();
}

//----------------------------------------------------------------------

void InputData::read_hierarchy 
(
 Hierarchy  * hierarchy,
 const FieldDescr * field_descr
 ) throw()
{

  IoHierarchy io_hierarchy(hierarchy);

  // Read hierarchy meta-data

  Input::read_meta (&io_hierarchy);

  // Call read_patch() on contained patches
  Input::read_hierarchy (hierarchy, field_descr);

}

//----------------------------------------------------------------------

Patch * InputData::read_patch 
(
 Patch * patch,
 const FieldDescr * field_descr,
 int ixp0, int iyp0, int izp0
 ) throw()
{
  ERROR("InputData::read_patch()",
	"Not Implemented");
  // if (patch == 0) {
  //   // create an uninitialized Patch
  //   patch = 
  //     factory_->create_patch
  //     (0,
  //      0,0,0,
  //      0,0,0,
  //      0,0,0,
  //      0,0,0,
  //      0,0,0,
  //      0);
  // }

  // Change to file group for patch

  char buffer[40];
  int id = patch->id();
  sprintf (buffer,"patch_%d",id);
  file_->group_chdir(buffer);
  file_->group_open();

  // Read patch meta-data

  IoPatch io_patch(patch);

  Input::read_meta_group (&io_patch);

  // Also read the patches parallel Layout

  IoLayout io_layout(patch->layout());

  Input::read_meta_group (&io_layout);

  // // Call read_block() on contained blocks
  // Input::read_patch(patch,field_descr,ixp0,iyp0,izp0);

  return patch;
}

//----------------------------------------------------------------------

void InputData::end_read_patch() throw()
{
  file_->group_close();
  file_->group_chdir("..");
}
//----------------------------------------------------------------------

Block * InputData::read_block 
( 
 Block * block,
 std::string  block_name,
 const FieldDescr * field_descr) throw()
{

  file_->group_chdir(block_name);
  file_->group_open();

  // Create temporary block

  ERROR("InputData::read_block",
	"Need to pass in Patch proxy to create_block");

  // block = factory_->create_block
  //   (0,0,0,
  //    0,0,0,
  //    0,0,0,
  //    0.0, 0.0, 0.0,
  //    0.0, 0.0, 0.0,
  //    1);

  // Read block meta data

  io_block()->set_block(block);

  Input::read_meta_group (io_block());

  int ibx=100,iby=100,ibz=100;
  block->index_patch(&ibx,&iby,&ibz);

  // // Call read_block() on base Input object

  // Input::read_block(block,field_descr,ixp0,iyp0,izp0);

  file_->group_close();
  file_->group_chdir("..");

  return block;
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

    // Read ith FieldBlock data

    file_->data_open(name.c_str(),&type,&nx,&ny,&nz);
    file_->data_read(buffer);
    file_->data_close();

    // Get ith FieldBlock data
    io_field_block()->data_value(i, &buffer, &name, &type, 
				 &nxd,&nyd,&nzd,
				 &nx, &ny, &nz);

  }

}

//======================================================================
