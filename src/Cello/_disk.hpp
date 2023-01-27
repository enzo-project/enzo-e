// See LICENSE_CELLO file for license and copyright information

/// @file     _disk.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-03-18 17:48:36
/// @brief    Private include file for the \ref Disk component

#ifndef _DISK_HPP
#define _DISK_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

// we explicitly avoid including <hdf5.h> here (including that header here
// forces it to be visible to all other dependent components)

#include <string>

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "disk_File.hpp"

// we explicitly avoid including "disk_FileHdf5.hpp" here and instead include
// in "disk_FileHdf5.cpp".
// -> Including that header here would force us to also include <hdf5.h> here &
//    make them both visible to all other dependent components
// -> if we decide to change this, then HDF5_C needs to be labelled as a public
//    linked library in CMakeLists.txt 

#endif /* _DISK_HPP */
