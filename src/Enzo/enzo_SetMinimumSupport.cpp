// See LICENSE_ENZO file for license and copyright information

/// @file      enzo_SetMinimumSupport.cpp
/// @author    Greg Bryan
/// @date      November, 1998
/// @brief     Set the energy to provide minimal pressure support

#include "cello.hpp"

#include "enzo.hpp"
 
//----------------------------------------------------------------------
 
int EnzoBlock::SetMinimumSupport(enzo_float &MinimumSupportEnergyCoefficient,
				 bool comoving_coordinates)
{
  if (NumberOfBaryonFields > 0) {
 
    const enzo_float pi = 3.14159;
 
    /* Compute cosmology factors. */
 
    enzo_float cosmo_a = 1.0, cosmo_dadt = 0.0;

    EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology * )
      simulation()->problem()->physics("cosmology");

    ASSERT ("EnzoBlock::SetMinimumSupport()",
	    "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
	    ! (comoving_coordinates && (cosmology == NULL)) );

    if (comoving_coordinates) {
      cosmology ->compute_expansion_factor(&cosmo_a, &cosmo_dadt,time());
    }

    enzo_float CosmoFactor = 1.0/cosmo_a;
 
    /* Determine the size of the grids. */
 
    const int in = cello::index_static();

    int dim, size = 1, i;
    for (dim = 0; dim < GridRank[in]; dim++)
      size *= GridDimension[dim];
 
    Field field = data()->field();

    enzo_float * density         = (enzo_float*) field.values("density");
    enzo_float * total_energy    = (enzo_float *)field.values("total_energy");
    enzo_float * internal_energy = (enzo_float *)field.values("internal_energy");
    enzo_float * velocity_x      = (enzo_float*) field.values("velocity_x");
    enzo_float * velocity_y      = (enzo_float*) field.values("velocity_y");
    enzo_float * velocity_z      = (enzo_float*) field.values("velocity_z");

    /* Set minimum GE. */
 
    MinimumSupportEnergyCoefficient =
      GravitationalConstant[in]/(4.0*pi) / (pi * (Gamma[in]*(Gamma[in]-1.0))) *
      CosmoFactor * MinimumPressureSupportParameter[in] *
      CellWidth[0] * CellWidth[0];
 
 
    /* PPM: set GE. */
 
    if (DualEnergyFormalism[in] == TRUE) {
      for (i = 0; i < size; i++)
	internal_energy[i] = MAX(internal_energy[i],
				 MinimumSupportEnergyCoefficient*density[i]);
      if (GridRank[in] != 3) return ENZO_FAIL;
      for (i = 0; i < size; i++)
	total_energy[i] = 
	  MAX((enzo_float)
	      (  internal_energy[i] +
		 0.5*(velocity_x[i]*velocity_x[i] +
		      velocity_y[i]*velocity_y[i] +
		      velocity_z[i]*velocity_z[i])),
	      total_energy[i]);
								
    }
    else {
      fprintf(stderr, "not implemented.\n");
      return ENZO_FAIL;
    }
 
  } // end: if (NumberOfBaryonFields > 0)
 
  return ENZO_SUCCESS;
}
