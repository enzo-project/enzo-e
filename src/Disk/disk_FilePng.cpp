// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FilePng.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @brief     Implementation of the FilePng class

#include "cello.hpp"

#include "disk.hpp"
 
//----------------------------------------------------------------------

FilePng::FilePng
(
 std::string path,
 std::string name
 ) throw() : File(path,name) 
{
  INCOMPLETE("FilePng::FilePng");
}

//----------------------------------------------------------------------
    
void FilePng::file_open () throw()
{
  INCOMPLETE("FilePng::file_open");
}

//----------------------------------------------------------------------

void FilePng::file_close () throw()
{
  INCOMPLETE("FilePng::file_close");
}

//----------------------------------------------------------------------

void FilePng::data_read 
( void * buffer, std::string name, enum scalar_type * type, 
  int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  INCOMPLETE("FilePng::data_read");
}

//----------------------------------------------------------------------

void FilePng::data_write 
( const void * buffer, std::string name, enum scalar_type type, 
  int n0, int n1, int n2, int n3, int n4) throw()
{
  INCOMPLETE("FilePng::data_read");
}

