// See LICENSE_CELLO file for license and copyright information

/// @file     disk_File.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-26
/// @brief    Implementation of the File class

#include "cello.hpp"

#include "disk.hpp"

//----------------------------------------------------------------------

File::File(std::string path, std::string name) throw ()
  : path_(path),
    name_(name)
{}
