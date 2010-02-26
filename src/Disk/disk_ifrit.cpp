/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @note      
 *
 *--------------------------------------------------------------------
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
 * @file      disk_ifrit.cpp
 * @brief     Implementation of the Ifrit class
 * @author    James Bordner
 * @date      Thu Feb 21 16:11:36 PST 2008
 *
 * $Id$
 *
 *********************************************************************
 */

#include <stdio.h>

#include "error.hpp"

#include "disk_ifrit.hpp"
 
//----------------------------------------------------------------------

void Ifrit::read_bin
(std::string name, Scalar * buffer, int * nx, int * ny, int * nz)
/**
 */
{
  FILE * fp = fopen (name.c_str(), "r");
  assert (sizeof(int)==4);
  assert (sizeof(float)==4);
  fread (nx,4,1,fp);
  fread (ny,4,1,fp);
  fread (nz,4,1,fp);
  fread (buffer,4,(*nx)*(*ny)*(*nz),fp);
}

//----------------------------------------------------------------------

void Ifrit::write_bin 
(std::string name, Scalar * buffer, int   nx, int   ny, int   nz)
/**
 */
{
  FILE * fp = fopen (name.c_str(), "w");
  assert (sizeof(int)==4);
  assert (sizeof(float)==4);
  fwrite (&nx,4,1,fp);
  fwrite (&ny,4,1,fp);
  fwrite (&nz,4,1,fp);
  fwrite (buffer,4,nx*ny*nz,fp);
}

