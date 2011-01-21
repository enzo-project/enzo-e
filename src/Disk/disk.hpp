// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef DISK_HPP
#define DISK_HPP

/// @file     disk.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-03-18 17:48:36
/// @brief    Include file for the \ref Disk component

#ifdef CONFIG_USE_HDF5
#  include <hdf5.h>
#endif

#include <string>

#include "cello.hpp"

#include "disk_File.hpp"
#include "disk_FileHdf5.hpp"
#include "disk_FileIfrit.hpp"

#endif /* DISK_HPP */
