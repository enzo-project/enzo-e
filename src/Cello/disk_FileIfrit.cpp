// See LICENSE_CELLO file for license and copyright information

/// @file      disk_FileIfrit.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @brief     Implementation of the FileIfrit class

#include "cello.hpp"

#include "disk.hpp"

//----------------------------------------------------------------------
 
FileIfrit::FileIfrit (std::string path, std::string name) throw()
  : File(path,name)
{
}

//----------------------------------------------------------------------

FileIfrit::~FileIfrit () throw()
{
}

//----------------------------------------------------------------------
 
void FileIfrit::file_open () throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::file_create () throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::file_close () throw()
{
}

//----------------------------------------------------------------------
  
void FileIfrit::file_read_meta
( void * buffer, std::string name,  int * s_type,
  int * n1, int * n2, int * n3, int * n4) throw()
{
}

//----------------------------------------------------------------------
  
void FileIfrit::file_write_meta
( const void * buffer, std::string name, int type,
  int n1, int n2, int n3, int n4) throw()
{
}

//----------------------------------------------------------------------
  

// Datasets

void FileIfrit::data_open
( std::string name,  int * type,
  int * m1, int * m2, int * m3, int * m4) throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::data_slice
( int   m1, int   m2, int   m3, int   m4,
  int   n1, int   n2, int   n3, int   n4,
  int   o1, int   o2, int   o3, int   o4) throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::data_create
( std::string name,  int type,
  int m1, int m2, int m3, int m4,
  int n1, int n2, int n3, int n4,
  int o1, int o2, int o3, int o4) throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::data_read (void * buffer) throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::data_write 
(const void * buffer) throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::data_close () throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::data_read_meta
( void * buffer, std::string name,  int * s_type,
  int * n1, int * n2, int * n3, int * n4) throw()
{
}

//----------------------------------------------------------------------
  
void FileIfrit::data_write_meta
( const void * buffer, std::string name, int type,
  int n1, int n2, int n3, int n4) throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::mem_create
  ( int mx, int my, int mz,
    int nx, int ny, int nz,
    int gx, int gy, int gz )
{
}

//----------------------------------------------------------------------

int FileIfrit::group_count () const throw()
{
  return 0;
}

//----------------------------------------------------------------------

std::string FileIfrit::group_name (size_t i) const throw()
{
  return "";
}

//----------------------------------------------------------------------

void FileIfrit::group_chdir (std::string name) throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::group_open () throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::group_create () throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::group_close () throw()
{
}

//----------------------------------------------------------------------

void FileIfrit::group_read_meta
( void * buffer, std::string name,  int * s_type,
  int * n1, int * n2, int * n3, int * n4) throw()
{
}

//----------------------------------------------------------------------
  
void FileIfrit::group_write_meta
( const void * buffer, std::string name, int type,
  int n1, int n2, int n3, int n4) throw()
{
}
