/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * SYNOPSIS:
 *
 *    
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * SYNOPSIS:
 *
 *    
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      image_dump.cpp
 * @brief     Project fields to images
 * @author    James Bordner
 * @date      Wed Nov 18 2009
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
 *    image_dump ();
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

#include "monitor.hpp"

void image_dump(const char * file_root, int cycle, double lower, double upper)
{ 

  int nx = GridDimension[0];
  int ny = GridDimension[1];
  int nz = GridDimension[2];

  Monitor monitor;

  // Open hdf5 file dump for cycle
  char filename[80];
  

  // color map
  double map[] = {1,1,1, 0,0,0};

  // slice
  sprintf (filename,"slice-%s-%06d-z.png",file_root,cycle);
  monitor.image(filename,
		BaryonField[field_density],nx,ny,nz,
		3,3,3,nx-3,ny-3,4,
		2,reduce_sum, lower/nx, upper/nx, map,2);
  
  if (nz > 1) {
    // projection
    sprintf (filename,"project-%s-%06d-z.png",file_root,cycle);
    monitor.image(filename,
		  BaryonField[field_density],nx,ny,nz,
		  3,3,3,nx-3,ny-3,nz-3,
		  2,reduce_sum,lower, upper, map,2);
  }

}
