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

  file_ = new FileHdf5 (".",file_name + ".h5");
  

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
( const Simulation * simulation ) throw()
/// Note factory, field_descr, and hierarchy needed since otherwise
/// Simulation functions must be called, which would introduce a circular
/// dependence between Simulation and Output components
{

  // Write parameter file

  bool is_root = simulation->group_process()->is_root();

  if (is_root) {

    std::string file_name = expand_file_name() + ".in";

    Parameters * parameters = simulation->parameters();

    // Update Initial parameters

    parameters->set_integer ("Initial:cycle",simulation->cycle());
    parameters->set_float   ("Initial:time", simulation->time());

    // Write restart parameter file
    parameters->write(file_name.c_str());
    
  }

  Output::write_simulation(simulation);

}

//----------------------------------------------------------------------
void OutputRestart::write_hierarchy 
(
 const Hierarchy * hierarchy,
 const FieldDescr * field_descr
 ) throw()
{
  Output::write_hierarchy(hierarchy,field_descr);
}

//----------------------------------------------------------------------

void OutputRestart::write_patch 
(
 const Patch * patch,
 const FieldDescr * field_descr,
 int ixp0, int iyp0, int izp0
 ) throw()
{
  Output::write_patch(patch,field_descr,ixp0,iyp0,izp0);
}

//----------------------------------------------------------------------

void OutputRestart::write_block 
(
 const Block * block,
 const FieldDescr * field_descr,
 int ixp0, int iyp0, int izp0) throw()
{
  Output::write_block(block,field_descr,ixp0,iyp0,izp0);
}

//----------------------------------------------------------------------

void OutputRestart::write_field
(
  const FieldBlock * field_block,
  const FieldDescr * field_descr,
  int field_index) throw()
{
  
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
