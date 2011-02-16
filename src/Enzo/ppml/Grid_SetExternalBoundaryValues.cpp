// $Id$
// See LICENSE_ENZO file for license and copyright information

/***********************************************************************
/
/  GRID CLASS (SET EXTERNAL BOUNDARY VALUES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    ENZO_SUCCESS or ENZO_FAIL
/
************************************************************************/
 
// Copy the current baryon fields to the old baryon fields/

#include "cello.hpp"

#include "enzo.hpp" 
 
int EnzoDescr::SetExternalBoundaryValues()
{
  int dim, field;
 
  /* Compute offset from corner of domain. */
 
  int GridOffset[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    if (dim < GridRank)
      GridOffset[dim] = NINT((GridLeftEdge[dim] - DomainLeftEdge[dim])/
			     CellWidth[dim]);
    else
      GridOffset[dim] = 0;
 
  /* loop through fields, setting each */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {
 
    if (BaryonField[field] == NULL) {
      fprintf(stderr, "Baryon field missing.\n");
      return ENZO_FAIL;
    }

#ifdef OOC_BOUNDARY
    ExternalBoundaryField = field;
#endif
 
    if (SetExternalBoundary(GridRank, GridDimension, GridOffset,
				      GridStartIndex, GridEndIndex,
				      BaryonField[field], FieldType[field])
	== ENZO_FAIL) {
      fprintf(stderr, "Error in SetExternalBoundary.\n");
      return ENZO_FAIL;
    }
 
  }
 
  return ENZO_SUCCESS;
 
}
