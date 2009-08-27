
#include "cello_hydro.h"

/* Cosmology Parameters */

/* The Hubble constant at z=0 (now) in units of 100 km/s/Mpc. */

float HubbleConstantNow;

/* The value of Omega due to non-relativistic particles at z=0. */

float OmegaMatterNow;

/* The value of Omega due to lamba (the cosmological constant) at z=0. */

float OmegaLambdaNow;

/* The comoving size of the simulation box (along the x-dir) in h^{-1} Mpc. */

float ComovingBoxSize;

/* The maximum allowed fractional increase of the expansion. */

float MaxExpansionRate;

/* The time in code units at the initial redshift. */

FLOAT InitialTimeInCodeUnits;

/* The initial redshift. */

FLOAT InitialRedshift;

/* Redshift output information: */

FLOAT CosmologyOutputRedshift[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
char *CosmologyOutputRedshiftName[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
FLOAT CosmologyOutputRedshiftTime[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
