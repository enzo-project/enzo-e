// See LICENSE_CELLO file for license and copyright information

/// @file     disk_File.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the File class

#include "cello.hpp"

#include "disk.hpp"

File::File(std::string path, std::string name, std::string mode) throw ()
  : path_(path),
    name_(name),
    mode_(mode)
{
}
