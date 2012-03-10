// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @brief    Implementation of the OutputData class

#include "cello.hpp"
#include "io.hpp"

//----------------------------------------------------------------------

OutputData::OutputData(const Factory * factory) throw ()
  : Output(factory)
{
}

//----------------------------------------------------------------------

OutputData::~OutputData() throw ()
{
}

//======================================================================

void OutputData::init () throw()
{
}

//----------------------------------------------------------------------

void OutputData::open () throw()
{
  std::string file_name = expand_file_name();

  Monitor::instance()->print ("Output","writing data file %s", 
			      file_name.c_str());

  delete file_;

  file_ = new FileHdf5 (".",file_name);

  file_->file_create();
}

//----------------------------------------------------------------------

void OutputData::close () throw()
{
  file_->file_close();

  delete file_;
  file_ = 0;
}

//----------------------------------------------------------------------

void OutputData::finalize () throw ()
{
  Output::finalize();
}

//----------------------------------------------------------------------

void OutputData::write_hierarchy 
(
 const Hierarchy * hierarchy,
 const FieldDescr * field_descr
 ) throw()
{

  // Loop over metadata items in Hierarchy

  IoHierarchy io_hierarchy(hierarchy);

  for (int i=0; i<io_hierarchy.meta_count(); i++) {

    std::string name;
    scalar_type type;
    void * buffer;
    int nx,ny,nz;

    // Get ith Hierarchy metadata
    io_hierarchy.meta_value(i,& buffer, &name, &type, &nx,&ny,&nz);

    // Write ith Hierarchy metadata
    file_->file_write_meta(buffer,name.c_str(),type,nx,ny,nz);
  }

  // Call write_patch() on contained patches
  Output::write_hierarchy (hierarchy, field_descr);

}

//----------------------------------------------------------------------

void OutputData::write_patch 
(
 const Patch * patch,
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

  for (int i=0; i<io_patch.meta_count(); i++) {

    void * buffer;
    std::string name;
    scalar_type type;
    int nx,ny,nz;

    // Get ith Patch metadata
    io_patch.meta_value(i,& buffer, &name, &type, &nx,&ny,&nz);

    // Write ith Patch metadata
    file_->group_write_meta(buffer,name.c_str(),type,nx,ny,nz);
  }

  // Call write_block() on contained blocks
  Output::write_patch(patch,field_descr,ixp0,iyp0,izp0);

  // BUG: this is getting executed before remote blocks begin writing
  // need to add a dependency so that this gets called only after last
  // block finishes

  file_->group_close();
  file_->group_chdir("..");

}

//----------------------------------------------------------------------

void OutputData::write_block 
( 
  const Block * block,
  const FieldDescr * field_descr,
  int ixp0, int iyp0, int izp0) throw()
{

  // Create file group for block

  char buffer[40];
  int ib = block->index();
  sprintf (buffer,"block_%d",ib);
  file_->group_chdir(buffer);
  file_->group_create();

  // Write block meta data

  io_block()->set_block(block);

  for (int i=0; i<io_block()->meta_count(); i++) {

    void * buffer;
    std::string name;
    scalar_type type;
    int nx,ny,nz;

    // Get ith Block metadata
    io_block()->meta_value(i,& buffer, &name, &type, &nx,&ny,&nz);

    // Write ith Block metadata
    file_->group_write_meta(buffer,name.c_str(),type,nx,ny,nz);

  }

  // Call write_block() on base Output object

  Output::write_block(block,field_descr,ixp0,iyp0,izp0);

  file_->group_close();
  file_->group_chdir("..");

}

//----------------------------------------------------------------------

void OutputData::write_field
( 
  const FieldBlock * field_block,
  const FieldDescr * field_descr,
  int field_index) throw()
{
  io_field_block()->set_field_descr(field_descr);
  io_field_block()->set_field_block(field_block);
  io_field_block()->set_field_index(field_index);

  for (int i=0; i<io_field_block()->data_count(); i++) {

    void * buffer;
    std::string name;
    scalar_type type;
    int nx,ny,nz;

    // Get ith FieldBlock data
    io_field_block()->data_value(i, &buffer, &name, &type, &nx,&ny,&nz);

    // Write ith FieldBlock data

    file_->data_create(name.c_str(),type,nx,ny,nz);
    file_->data_write(buffer);
    file_->data_close();
  }

}

//----------------------------------------------------------------------

void OutputData::prepare_remote (int * n, char ** buffer) throw()
{
}

//----------------------------------------------------------------------

void OutputData::update_remote  ( int n, char * buffer) throw()
{
}

//----------------------------------------------------------------------


void OutputData::cleanup_remote (int * n, char ** buffer) throw()
{
}

//======================================================================
