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
{
}

//----------------------------------------------------------------------

void File::file_read_meta
( void * buffer, std::string name,  enum scalar_type * s_type,
  int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  WARNING("File::file_read_meta",
	  "Operation not supported by object");
}

//----------------------------------------------------------------------
  
void File::file_write_meta
( const void * buffer, std::string name, enum scalar_type type,
  int n0, int n1, int n2, int n3, int n4) throw()
{
  WARNING("File::file_write_meta",
	  "Operation not supported by object");
}

//----------------------------------------------------------------------

void File::data_read_meta
( void * buffer, std::string name,  enum scalar_type * s_type,
  int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  WARNING("File::data_read_meta",
	  "Operation not supported by object");
}

//----------------------------------------------------------------------
  
void File::data_write_meta
( const void * buffer, std::string name, enum scalar_type type,
  int n0, int n1, int n2, int n3, int n4) throw()
{
  WARNING("File::data_write_meta",
	  "Operation not supported by object");
}

//----------------------------------------------------------------------

void File::group_open (std::string name) throw()
{
  WARNING("File::group_open",
	  "Operation not supported by object");
}

//----------------------------------------------------------------------

void File::group_create (std::string name) throw()
{
  WARNING("File::group_create",
	  "Operation not supported by object");
}

//----------------------------------------------------------------------

void File::group_close () throw()
{
  WARNING("File::group_close",
	  "Operation not supported by object");
}

//----------------------------------------------------------------------

void File::group_read_meta
( void * buffer, std::string name,  enum scalar_type * s_type,
  int * n0, int * n1, int * n2, int * n3, int * n4) throw()
{
  WARNING("File::group_read_meta",
	  "Operation not supported by object");
}

//----------------------------------------------------------------------
  
void File::group_write_meta
( const void * buffer, std::string name, enum scalar_type type,
  int n0, int n1, int n2, int n3, int n4) throw()
{
  WARNING("File::group_write_meta",
	  "Operation not supported by object");
}

