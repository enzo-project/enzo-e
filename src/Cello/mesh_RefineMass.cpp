// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineMass.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-23
/// @brief    Implementation of RefineMass class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineMass::RefineMass
(
 double min,
 double level_exponent,
 double min_overdensity,
 double root_cell_volume) throw ()
  : min_(min),
    level_exponent_(level_exponent),
    min_overdensity_(min_overdensity)
  // ENZO non-cosmology

  //    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
  //    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
  //      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
  //        float(MetaData.TopGridDims[dim]);

  // ENZO Cosmology
  //      MinimumMassForRefinement[i] = CosmologySimulationOmegaBaryonNow/
  //	                            OmegaMatterNow;
  //      if (CellFlaggingMethod[i] == 4)
  //	MinimumMassForRefinement[i] = CosmologySimulationOmegaCDMNow/
  //	                              OmegaMatterNow;
  // 
  //      MinimumMassForRefinement[i] *= MinimumOverDensityForRefinement[i];
  //      for (dim = 0; dim < MetaData.TopGridRank; dim++)
  //	MinimumMassForRefinement[i] *=
  //	  (DomainRightEdge[dim]-DomainLeftEdge[dim])/
  //	  float(MetaData.TopGridDims[dim]);
{
  TRACE("RefineSlope::RefineSlope");
  WARNING ("RefineMass::RefineMass()",
	   "Assuming non-Cosmology problem for RefineMass");
  if (min_ == -1.0) {
    min_ = min_overdensity * root_cell_volume;
  }

}

//----------------------------------------------------------------------

int RefineMass::apply 
(
 CommBlock * comm_block,
 const FieldDescr * field_descr
 ) throw ()
{
  FieldBlock * field_block = comm_block->block()->field_block();
  int level = 0;

  WARNING ("RefineMass::RefineMass()",
	   "Assuming refinement factor == 2");
  WARNING ("RefineMass::RefineMass()",
	   "Assuming level==0 for level_exponent");

  double mass_min_refine  = min_ * pow(2.0,level*level_exponent_);
  double mass_max_coarsen = mass_min_refine / 2.0;

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  int gx,gy,gz;
  field_descr->ghosts(0, &gx,&gy,&gz);

  precision_type precision = field_descr->precision(0);

  void * void_array = field_block->field_values(0);

  //  int num_fields = field_descr->field_count();

  bool all_coarsen = true;
  bool any_refine = false;

  const int d3[3] = {1,nx,nx*ny};

  switch (precision) {
  case precision_single:
    {
      float * array = (float*)void_array;
      float mass;
      for (int axis=0; axis<3; axis++) {
	int d = d3[axis];
	for (int ix=0; ix<nx; ix++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int iz=0; iz<nz; iz++) {
	      int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
	      mass = (array[i]) ? (array[i+d] - array[i-d]) / array[i] : 0.0;
	      if (mass > mass_min_refine)  any_refine  = true;
	      if (mass > mass_max_coarsen) all_coarsen = false;
	    }
	  }
	}
      }
    }
    break;
  case precision_double:
    {
      double * array = (double*)void_array;
      double mass;

      for (int axis=0; axis<3; axis++) {
	int d = d3[axis];
	for (int ix=0; ix<nx; ix++) {
	  for (int iy=0; iy<ny; iy++) {
	    for (int iz=0; iz<nz; iz++) {
	      int i = (gx+ix) + nx*((gy+iy) + ny*(gz+iz));
	      mass = (array[i]) ? (array[i+d] - array[i-d]) / array[i] : 0.0;
	      if (mass > mass_min_refine)  any_refine  = true;
	      if (mass > mass_max_coarsen) all_coarsen = false;
	    }
	  }
	}
      }
    }
    break;
  default:
    ERROR2("RefineMass::apply",
	   "Unknown precision %d for field %d",
	   precision,0);
    break;
  }

  return 
    any_refine ?  adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

}


//======================================================================

