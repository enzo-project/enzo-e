// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @todo     Move write_hierarchy() call to write_patch() etc. to parent
///           Output::write_hierarchy()
/// @brief    Implementation of the OutputData class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

OutputData::OutputData(Simulation * simulation) throw ()
  : Output(simulation)
{
}

//----------------------------------------------------------------------

OutputData::~OutputData() throw ()
{
}

//======================================================================

void OutputData::init () throw()
{
  TRACE("OutputData::init ()");
}

//----------------------------------------------------------------------

void OutputData::open () throw()
{
  TRACE("OutputData::open()");

  std::string file_name = expand_file_name();

  Monitor::instance()->print ("[Output] writing data file %s", 
			      file_name.c_str());

  delete file_;

  file_ = new FileHdf5 (".",file_name);

  file_->file_create();
}

//----------------------------------------------------------------------

void OutputData::close () throw()
{
  TRACE("OutputData::close ()");
  file_->file_close();

  delete file_;
  file_ = 0;
}

//----------------------------------------------------------------------

void OutputData::finalize () throw ()
{
  TRACE("OutputData::finalize ()");
  Output::finalize();
}

//----------------------------------------------------------------------

void OutputData::write_hierarchy 
(
 const FieldDescr * field_descr,
 Hierarchy * hierarchy
 ) throw()
{
  TRACE("ENTER OutputData::write_hierarchy ()");

  // (*)  write meta factory_

  /* initialized in constructor */

  // (*)  write meta patch_list_

  int patch_count = hierarchy->num_patches();
  file_->file_write_meta(&patch_count,"patch_count",scalar_type_int);
  
  // (*)  write meta tree_

  /* @@@ */

  // (*)  write meta lower_[];

  double lower[3];
  hierarchy->lower(&lower[0],&lower[1],&lower[2]);
  file_->file_write_meta(lower,"domain_lower",scalar_type_double,3);

  // (*)  double meta upper_[3];

  double upper[3];
  hierarchy->upper(&upper[0],&upper[1],&upper[2]);
  file_->file_write_meta(upper,"domain_upper",scalar_type_double,3);

  ItPatch it_patch (hierarchy);

  // (*) write data patch_list_

  while (Patch * patch = ++it_patch) {

    // NO OFFSET: ASSUMES ROOT PATCH
    write_patch (field_descr, patch,  0,0,0);

  }

  // ( ) write data meta tree_
  TRACE("EXIT OutputData::write_hierarchy ()");

}

//----------------------------------------------------------------------

void OutputData::write_patch 
(
 const FieldDescr * field_descr,
 Patch * patch,
 int ixp0, int iyp0, int izp0
 ) throw()
{
  TRACE("ENTER OutputData::write_patch ()");

  INCOMPLETE("OutputData::write_patch(): write patch metadata");

  Output::write_patch(field_descr,patch,ixp0,iyp0,izp0);
  
  TRACE("EXIT OutputData::write_patch ()");

}

//----------------------------------------------------------------------

void OutputData::write_block ( const FieldDescr * field_descr,
  Block * block,
  int ixp0, int iyp0, int izp0) throw()
{
  // Get block index

  int ib = block->index();

  char buffer[40];
  sprintf (buffer,"/block-%d",ib);

  TRACE1 ("OutputData::write_block() %s",buffer);
  file_->group_create(buffer);

  file_->group_write_meta(&ib,"block_index",scalar_type_int);

  TRACE("ENTER OutputData::write_block ()");

  INCOMPLETE("OutputData::write_block(): write block metadata");
  INCOMPLETE("OutputData::write_block(): write block fields");
  TRACE("EXIT OutputData::write_block ()");
  file_->group_close();

}

//----------------------------------------------------------------------

void OutputData::prepare_remote (int * n, char ** buffer) throw()
{
  TRACE("OutputData::prepare_remote ()");
}

//----------------------------------------------------------------------

void OutputData::update_remote  ( int n, char * buffer) throw()
{
  TRACE("OutputData::update_remote  ()");
}

//----------------------------------------------------------------------


void OutputData::cleanup_remote (int * n, char ** buffer) throw()
{
  TRACE("OutputData::cleanup_remote ()");
}

//======================================================================
