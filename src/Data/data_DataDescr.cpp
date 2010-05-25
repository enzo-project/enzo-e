// $Id: data_DataDescr.cpp 1388 2010-04-20 23:57:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     data_DataDescr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the DataDescr class

#include "data_DataDescr.hpp"

DataDescr::DataDescr() throw ()
  : field_descr_(),
    particle_descr_()
{ 
}

//----------------------------------------------------------------------

FieldDescr * DataDescr::field_descr ()
{ 
  return & field_descr_;
}

