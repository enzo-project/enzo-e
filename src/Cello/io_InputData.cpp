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

//----------------------------------------------------------------------

InputData::~InputData() throw()
{
  close();
}

//----------------------------------------------------------------------

void InputData::pup (PUP::er &p)
{

  TRACEPUP;

  // NOTE: change this function whenever attributes change

  Input::pup(p);

  // this function intentionally left blank
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

  // Calls read_blocks() on contained forest
  Input::read_hierarchy (hierarchy, field_descr);

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

  // Read block meta data

  io_block()->set_block(block);

  Input::read_meta_group (io_block());

  int ibx,iby,ibz;
  block->index_forest(&ibx,&iby,&ibz);

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
  io_field_block()->set_field_descr((FieldDescr*)field_descr);
  io_field_block()->set_field_block(field_block);
  io_field_block()->set_field_index(field_index);

  for (size_t i=0; i<io_field_block()->data_count(); i++) {

    void * buffer = 0;
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
