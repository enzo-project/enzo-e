// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      data_dump.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Sun Aug 30 14:16:29 PDT 2009
/// @brief     Write BaryonField's to disk

#include "cello.hpp"

#include "enzo.hpp"

void data_dump(const char * file_root, int cycle)
{ 

  int nx = GridDimension[0];
  int ny = GridDimension[1];
  int nz = GridDimension[2];

  FileHdf5 file;

  // Open hdf5 file dump for cycle
  char filename[80];

  sprintf (filename,"%s-%06d",file_root,cycle);
  file.file_open(filename,"w");
  file.dataset_open_write ("density",nx,ny,nz);
  file.write(BaryonField[field_density]);
  file.dataset_close ();
  file.file_close();
}
