// See LICENSE_ENZO file for license and copyright information

/* foflib.h
   Mark Krumholz, 3/20/00
   Modified by Mark Krumholz, 8/22/00
   Modfiied by Nathan Goldbaum, December 2011, to include in enzo
   Header file to accompany foflib.c */

//#include "ErrorExceptions.h"
//#include "macros_and_parameters.h"
//#include "typedefs.h"
//#include "global_data.h"

int FofVar(int, enzo_float *, enzo_float *, int *, int **);
int FofVarList(int, enzo_float *, enzo_float *, int *, int **, int ***);
int Fof(int, enzo_float *, enzo_float, int *, int **);
int FofList(int, enzo_float *, enzo_float, int *, int **, int ***);
int FofPrune(int, int, int *, int **, int);
int FofListPrune(int, int, int *, int **, int ***, int);
void NearNeighbor(int, enzo_float *, int, int *);
void NearNeighborPartial(int, enzo_float *, int, int, int *, int *);
void FindNeighbor(int, enzo_float *, enzo_float, int ***, int *);
void FindNeighborPartial(int, enzo_float *, int, int *, enzo_float *, int ***, int *);
