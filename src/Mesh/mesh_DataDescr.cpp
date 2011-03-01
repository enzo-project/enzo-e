// $Id: template.cpp 2035 2011-02-28 23:47:31Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_DataDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file mesh_DataDescr.cpp
///
/// Detailed description of file mesh_DataDescr.cpp

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

DataDescr::DataDescr() throw ()
  : field_descr_(new FieldDescr)
{

}

//----------------------------------------------------------------------

DataDescr::~DataDescr() throw ()
{
  delete field_descr_;
}

//----------------------------------------------------------------------

DataDescr::DataDescr(const DataDescr & data_descr) throw ()
/// @param     DataDescr  Object being copied
  : field_descr_(new FieldDescr(*((const FieldDescr *)(data_descr.field_descr()))))
{
}

//----------------------------------------------------------------------

DataDescr & DataDescr::operator= (const DataDescr & data_descr) throw ()
/// @param     DataDescr  Source object of the assignment
/// @return    The target assigned object
{
  field_descr_ = new FieldDescr;
  (*field_descr_) = *(data_descr.field_descr());
  return *this;
}

//======================================================================





