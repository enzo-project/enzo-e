// See LICENSE_CELLO file for license and copyright information

/// @file      print_field.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Sat Aug 29 23:45:03 PDT 2009
/// @brief     Print the active values of the grid field

#include "cello.hpp"

#include "enzo.hpp"
 
void EnzoBlock::print_field (int field)

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
    
