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

#ifdef CONFIG_USE_CHARM

void OutputData::init 
(
 const Hierarchy * hierarchy,
 int cycle,
 double time
 ) throw()
{
  INCOMPLETE("OutputData::init()");
}

//----------------------------------------------------------------------

void OutputData::block (const Block * block) throw()
{
  INCOMPLETE("OutputData::block()");
}

#endif

//----------------------------------------------------------------------

void OutputData::write
(
 const FieldDescr * field_descr,
 Hierarchy * hierarchy,
 int cycle,
 double time,
 bool root_call
  ) throw()
{

  if (root_call) {
  }

  ItPatch it_patch (hierarchy);
  while (Patch * patch = ++it_patch) {
    write (field_descr, patch, cycle,time, false,  0,0,0);
  }

  if (root_call) {
  }

}

//----------------------------------------------------------------------

void OutputData::write
(
 const FieldDescr * field_descr,
 Patch * patch,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
 ) throw()
{
  if (root_call) {
  }

  ItBlock it_block (patch);
  while (Block * block = ++it_block) {
    write (field_descr, block, cycle,time, false,  0,0,0);
  }

  if (root_call) {
  }

}

//----------------------------------------------------------------------

void OutputData::write
(
 const FieldDescr * field_descr,
 Block * block,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
) throw()
{
  if (root_call) {
  }

  if (root_call){
  }

}

//======================================================================
