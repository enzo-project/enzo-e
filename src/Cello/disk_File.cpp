// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FileHdf5.cpp
/// @author    Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date      Jan 13 2023
/// @brief     Implementation of the File class's static functions/attributes

#include "cello.hpp"

#include "disk.hpp"

// we intentionally include "disk_FileHdf5.hpp" here, rather than in "disk.hpp",
// "_disk.hpp" or any other header included in the files so that dependents of
// Cello's disk component don't have a transitive dependence on the hdf5 library
#include "disk_FileHdf5.hpp"

std::map<const std::string,File *> File::file_list;

//----------------------------------------------------------------------

File* File::construct_FileHdf5(std::string path, std::string name) throw()
{ return new FileHdf5(path, name); }
