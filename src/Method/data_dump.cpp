//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */


/** 
 *********************************************************************
 *
 * @file      data_dump.cpp
 * @brief     Write BaryonField's to disk
 * @author    James Bordner
 * @date      Sun Aug 30 14:16:29 PDT 2009
 *
 * DESCRIPTION 
 * 
 *    Write BaryonField's to disk
 *
 * PACKAGES
 *
 *    NONE
 * 
 * INCLUDES
 *  
 *    cello_hydro.h
 *
 * PUBLIC FUNCTIONS
 *  
 *    data_dump ();
 *
 * PRIVATE FUCTIONS
 *  
 *    NONE
 *
 * $Id$
 *
 *********************************************************************
 */

#include "cello_hydro.h"

#include <string>
#include <hdf5.h>

#include "disk.hpp"

void data_dump(const char * file_root, int cycle)
{ 

  int nx = GridDimension[0];
  int ny = GridDimension[1];
  int nz = GridDimension[2];

  Hdf5 hdf5;

  // Open hdf5 file dump for cycle
  char filename[80];

  sprintf (filename,"%s-%06d.hdf5",file_root,cycle);
  hdf5.file_open(filename,"w");
  hdf5.dataset_open_write ("density",nx,ny,nz);
  hdf5.write(BaryonField[field_density]);
  hdf5.dataset_close ();
  hdf5.file_close();
}
