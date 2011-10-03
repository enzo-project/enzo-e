// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
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
  std::string file_name = expand_file_name();
  file_ = new FileHdf5 (".",file_name);
}

//----------------------------------------------------------------------

void OutputData::open () throw()
{
  TRACE("OutputData::open()");
  std::string file_name = expand_file_name();

  Monitor::instance()->print ("[Output] writing data file %s", 
			      file_name.c_str());

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
  TRACE("OutputData::write_hierarchy ()");
  TRACE("()");

  ItPatch it_patch (hierarchy);
  while (Patch * patch = ++it_patch) {
    write_patch (field_descr, patch,  0,0,0);
  }

}

//----------------------------------------------------------------------

void OutputData::write_patch 
(
 const FieldDescr * field_descr,
 Patch * patch,
 int ixp0, int iyp0, int izp0
 ) throw()
{
  TRACE("OutputData::write_patch ()");

  ItBlock it_block (patch);
  while (Block * block = ++it_block) {
    write_block (field_descr, block, 0,0,0);
  }

}

//----------------------------------------------------------------------

void OutputData::write_block ( const FieldDescr * field_descr,
  Block * block,
  int ixp0, int iyp0, int izp0) throw()
{
  TRACE("OutputData::write_block ()");

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
