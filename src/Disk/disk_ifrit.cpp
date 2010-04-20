// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      disk_ifrit.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 16:11:36 PST 2008
/// @brief     Implementation of the Ifrit class

#include <stdio.h>

#include "cello.h"
#include "error.hpp"
#include "disk_ifrit.hpp"
 
//----------------------------------------------------------------------

void Ifrit::read_bin
(std::string name, 
 Scalar *    buffer, 
 int *       nx, 
 int *       ny, 
 int *       nz) throw()
///
{
  FILE * fp = fopen (name.c_str(), "r");
  ASSERT ("Ifrit::read_bin","int or float size in unexpected",
	  sizeof(int)==4 && sizeof(float)==4);
  fread (nx,4,1,fp);
  fread (ny,4,1,fp);
  fread (nz,4,1,fp);
  fread (buffer,4,(*nx)*(*ny)*(*nz),fp);
}

//----------------------------------------------------------------------

void Ifrit::write_bin 
(std::string name, 
 Scalar *    buffer, 
 int         nx, 
 int         ny, 
 int         nz) throw()
///
{
  FILE * fp = fopen (name.c_str(), "w");
  ASSERT ("Ifrit::read_bin","int or float size in unexpected",
	  sizeof(int)==4 && sizeof(float)==4);
  fwrite (&nx,4,1,fp);
  fwrite (&ny,4,1,fp);
  fwrite (&nz,4,1,fp);
  fwrite (buffer,4,nx*ny*nz,fp);
}

