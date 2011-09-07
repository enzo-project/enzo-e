// See LICENSE_CELLO file for license and copyright information

/// @file     io_Io.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-06
/// @brief    Implementation of Io abstract base class functions
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

Io::Io(std::string directory,
   std::string filename,
   enum file_format_type file_format
   ) throw()
  : file_(NULL)
{

  is_file_new_ = true;

  switch (file_format) {

  case file_format_hdf5:

    file_ = new FileHdf5 (directory,filename);
    break;

  default:

    char message[ERROR_LENGTH];
    sprintf (message,"Unsupported file_format_type %d", file_format);
    WARNING("Io::Io",message);
    break;

  }
}

//----------------------------------------------------------------------

Io::Io(File & file) throw()
  : file_(&file)
{
  is_file_new_ = false;
}

//----------------------------------------------------------------------

Io::~Io() throw ()
{
  if (is_file_new_) delete file_;
}

//======================================================================

