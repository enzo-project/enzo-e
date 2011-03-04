// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_HPP
#define DISK_HPP

/// @file     disk.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-03-18 17:48:36
/// @todo     Remove CONFIG_USE_HDF5
/// @brief    Include file for the \ref Disk component

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>

#ifdef CONFIG_USE_HDF5
#  include <hdf5.h>
#endif

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "error.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "disk_File.hpp"
#include "disk_FileHdf5.hpp"
#include "disk_FileIfrit.hpp"

#endif /* DISK_HPP */
