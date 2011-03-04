// $Id$
// See LICENSE_ENZO file for license and copyright information

/***********************************************************************
/
/  GRID CLASS (IDENTIFY CERTAIN COMMONLY USED VARIABLES FROM THE LIST)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/

#include "cello.hpp"

#include "enzo.hpp"
 
/* function prototypes */
 
// int FindField(int f, int farray[], int n);
 
 
 
int EnzoBlock::IdentifyPhysicalQuantities
(
 int &DensNum, int &GENum, int &Vel1Num,
 int &Vel2Num, int &Vel3Num, int &TENum )
{
 
  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = 0;
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "GIPQ: Cannot find density.\n");
    return ENZO_FAIL;
  }
 
  /* Find Total energy, if possible. */
 
  if ((TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find total energy.\n");
    return ENZO_FAIL;
  }
 
  /* Find gas energy, if possible. */
 
  if (DualEnergyFormalism == TRUE)
    if ((GENum = FindField(InternalEnergy, FieldType,
			   NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find gas energy.\n");
      return ENZO_FAIL;
    }
 
  /* Find Velocity1, if possible. */
 
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find Velocity1.\n");
    return ENZO_FAIL;
  }
 
  /* Find Velocity2, if possible. */
 
  if (GridRank > 1)
    if ((Vel2Num = FindField(Velocity2, FieldType,
			     NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Cannot find Velocity2.\n");
      return ENZO_FAIL;
    }
 
  /* Find Velocity3, if possible. */
 
  if (GridRank > 2)
    if ((Vel3Num = FindField(Velocity3, FieldType,
			     NumberOfBaryonFields)) == 0) {
      fprintf(stderr, "Cannot find Velocity3.\n");
      return ENZO_FAIL;
    }
 
  return ENZO_SUCCESS;
}
