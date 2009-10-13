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

#include "array.hpp"
#include "disk.hpp"

void data_dump(const char * file_root, int cycle)
{ 

  Hdf5 hdf5;

  // Open hdf5 file dump for cycle
  char filename[80];
  sprintf (filename,"%s-%06d.hdf5",file_root,cycle);
  hdf5.file_open(filename,"w");

  // Open density dataset
  Array density (BaryonField[field_density],
		 GridDimension[0],
		 GridDimension[1],
		 GridDimension[2]);
  hdf5.dataset_open ("density",density);

  // Write the data to disk
  hdf5.write(density);

  // Close the dateset
  hdf5.dataset_close ();

  // Close the file
  hdf5.file_close();
}
