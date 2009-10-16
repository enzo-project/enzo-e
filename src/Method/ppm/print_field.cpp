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
 * @file      print_field.cpp
 * @brief     Print the active values of the grid field
 * @author    James Bordner
 * @date      Sat Aug 29 23:45:03 PDT 2009
 *
 * DESCRIPTION 
 * 
 *    Print the active values of the grid field
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
 *    NONE
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
 
void print_field (int field)

{
  if (BaryonField[field] == NULL) return -1;

  int ndx = GridDimension[0];
  int ndy = GridDimension[1];

  char buffer [80];
  sprintf (buffer,"cello-field-%d",field);
  FILE *fp = fopen (buffer,"w");

  for (int iz = GridStartIndex[2]; iz<=GridEndIndex[2]; iz++) {
    for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {
      for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {
	int i = ix + ndx * (iy + ndy * iz);
	fprintf (fp,"%d %d %d %g\n",ix,iy,iz,BaryonField[field][i]);
      }
    }
  }
  fclose(fp);

}
    
