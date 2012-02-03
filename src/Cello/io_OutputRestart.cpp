// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputRestart.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-02-02
/// @brief    Implementation of writing restart dumps
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

OutputRestart::OutputRestart(const Factory * factory) throw ()
  : Output(factory)
{
}

//----------------------------------------------------------------------

OutputRestart::~OutputRestart() throw ()
{
}

//======================================================================

void OutputRestart::init () throw()
{
}

//----------------------------------------------------------------------

void OutputRestart::open () throw()
{
  std::string file_name = expand_file_name();

  Monitor::instance()->print ("Output","writing data file %s", 
			      file_name.c_str());

  delete file_;

  file_ = new FileHdf5 (".",file_name);

  file_->file_create();
}

//----------------------------------------------------------------------

void OutputRestart::close () throw()
{
  file_->file_close();

  delete file_;
  file_ = 0;
}

//----------------------------------------------------------------------

void OutputRestart::finalize () throw ()
{
  Output::finalize();
}

//----------------------------------------------------------------------

void OutputRestart::write_simulation
(
 Factory    * factory,
 FieldDescr * field_descr,
 Hierarchy  * hierarchy,
 Simulation * simulation
 ) throw()
/// Note factory, field_descr, and hierarchy needed since otherwise
/// Simulation functions must be called, which would introduce a circular
/// dependence between Simulation and Output components
{
  OutputData output_data (factory);

  output_data.write_hierarchy(field_descr,hierarchy);
  
}

//----------------------------------------------------------------------
void OutputRestart::write_hierarchy 
(
 const FieldDescr * field_descr,
 Hierarchy * hierarchy
 ) throw()
{
  ERROR("OutputRestart::write_field()",
	"This function is not supported for this output type");
}

//----------------------------------------------------------------------

void OutputRestart::write_patch 
(
 const FieldDescr * field_descr,
 Patch * patch,
 int ixp0, int iyp0, int izp0
 ) throw()
{
  ERROR("OutputRestart::write_field()",
	"This function is not supported for this output type");
}

//----------------------------------------------------------------------

void OutputRestart::write_block ( const FieldDescr * field_descr,
  Block * block,
  int ixp0, int iyp0, int izp0) throw()
{
  ERROR("OutputRestart::write_field()",
	"This function is not supported for this output type");
}

//----------------------------------------------------------------------

void OutputRestart::write_field
( const FieldDescr * field_descr,
  FieldBlock * field_block,
  int field_index) throw()
{
  ERROR("OutputRestart::write_field()",
	"This function is not supported for this output type");
}

//----------------------------------------------------------------------

void OutputRestart::prepare_remote (int * n, char ** buffer) throw()
{
}

//----------------------------------------------------------------------

void OutputRestart::update_remote  ( int n, char * buffer) throw()
{
}

//----------------------------------------------------------------------


void OutputRestart::cleanup_remote (int * n, char ** buffer) throw()
{
}

//======================================================================
