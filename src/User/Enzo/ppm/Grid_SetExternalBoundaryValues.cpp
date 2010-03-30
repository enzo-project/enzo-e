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
/    SUCCESS or FAIL
/
************************************************************************/
 
// Copy the current baryon fields to the old baryon fields/
 
#include "cello_hydro.h"
 
int SetExternalBoundaryValues()
{
  int dim, field;
 
  /* Compute offset from corner of domain. */
 
  int GridOffset[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    if (dim < GridRank)
      GridOffset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/
			     CellWidth[dim][0]);
    else
      GridOffset[dim] = 0;
 
  /* loop through fields, setting each */
 
  for (field = 0; field < NumberOfBaryonFields; field++) {
 
    if (BaryonField[field] == NULL) {
      fprintf(stderr, "Baryon field missing.\n");
      return FAIL;
    }

#ifdef OOC_BOUNDARY
    ExternalBoundaryField = field;
#endif
 
    if (SetExternalBoundary(GridRank, GridDimension, GridOffset,
				      GridStartIndex, GridEndIndex,
				      BaryonField[field], FieldType[field])
	== FAIL) {
      fprintf(stderr, "Error in SetExternalBoundary.\n");
      return FAIL;
    }
 
  }
 
  return SUCCESS;
 
}
