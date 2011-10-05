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
 Hierarchy * hierarchy,
 int index_output_charm
 ) throw()
{
  TRACE("OutputData::write_hierarchy ()");

  INCOMPLETE("OutputData::write_hierarchy(): write hierarchy metadata");

  // ( )  const Factory * factory_;
  // ( )  std::vector<Patch *> patch_list_;
  // ( )  TreeK * tree_;

  // (*)  double lower_[3];
  double lower[3];
  hierarchy->lower(&lower[0],&lower[1],&lower[2]);
  file_->file_write_meta(lower,"lower",scalar_type_double,3);

  // (*)  double upper_[3];
  double upper[3];
  hierarchy->upper(&upper[0],&upper[1],&upper[2]);
  file_->file_write_meta(upper,"upper",scalar_type_double,3);

  ItPatch it_patch (hierarchy);

  while (Patch * patch = ++it_patch) {

      if (patch->blocks_allocated()) {
#ifdef CONFIG_USE_CHARM
	patch->block_array().p_write (index_output_charm);
#else
	// NO OFFSET: ASSUMES ROOT PATCH
	write_patch (field_descr, patch,  0,0,0);
#endif
      }
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

  INCOMPLETE("OutputData::write_patch(): write patch metadata");

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

  INCOMPLETE("OutputData::write_block(): write block metadata");
  INCOMPLETE("OutputData::write_block(): write block fields");

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
