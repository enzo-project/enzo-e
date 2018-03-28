// // See LICENSE_CELLO file for license and copyright information

// /// @file     enzo_EnzoRefineMass.cpp
// /// @author   James Bordner (jobordner@ucsd.edu)
// /// @date     2013-04-23
// /// @brief    Implementation of Enzo RefineMass class

// #define DEBUG_REFINE

#include "enzo.hpp"
#include "charm_simulation.hpp"
#include "enzo.hpp"
#include "enzo.decl.h"

//----------------------------------------------------------------------

EnzoRefineMass::EnzoRefineMass
(
 double min_refine,
 double max_coarsen,
 int    max_level,
 bool include_ghosts,
 std::string output,
 std::string name,
 std::string mass_type,
 double level_exponent) throw ()
  : Refine(min_refine,max_coarsen,max_level,include_ghosts,output),
    name_(name),
    mass_ratio_(0.0),
    level_exponent_(level_exponent)

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
  EnzoSimulation * simulation =
    (EnzoSimulation * ) proxy_simulation.ckLocalBranch();

  EnzoPhysicsCosmology * cosmology =
    (EnzoPhysicsCosmology *) simulation->problem()->physics("cosmology");

  if (cosmology) {
    if (mass_type == "dark") {
      mass_ratio_ = cosmology->omega_cdm_now()   /cosmology->omega_matter_now();
    } else if (mass_type == "baryon") {
      mass_ratio_ = cosmology->omega_baryon_now()/cosmology->omega_matter_now();
    } else {
      ERROR1 ("EnzoRefineMass::EnzoRefineMass()",
	      "Unknown mass_type %s",  mass_type.c_str());
    }
  } else {
    mass_ratio_ = 0.0;
  }
}

//----------------------------------------------------------------------

int EnzoRefineMass::apply ( Block * block ) throw ()
{
  Field field = block->data()->field();
  int level = block->level();

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  EnzoSimulation * simulation =
    (EnzoSimulation * ) proxy_simulation.ckLocalBranch();
  EnzoPhysicsCosmology * cosmology =
    (EnzoPhysicsCosmology *) simulation->problem()->physics("cosmology");

  double hx0 = hx*pow(2.0,level);
  double hy0 = hy*pow(2.0,level);
  double hz0 = hz*pow(2.0,level);
  double scale = (mass_ratio_ == 0.0) ? 1.0 :
    mass_ratio_*pow(2.0,level*level_exponent_)*hx0*hy0*hz0;
  double mass_min_refine  = scale*min_refine_;
  double mass_max_coarsen = scale*max_coarsen_;

  const int id_field = field.field_id(name_);
  ASSERT1 ("EnzoRefineMass::apply()",
	   "Undefined field name %s",
	   name_.c_str(), id_field >= 0);
  
  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (id_field, &mx,&my,&mz);
  field.ghost_depth(id_field, &gx,&gy,&gz);

  precision_type precision = field.precision(id_field);

  //  int num_fields = field_descr->field_count();

  bool all_coarsen = true;
  bool any_refine = false;

  union {
    float *       rho4;
    double *      rho8;
    long double * rho16;
  };
  union {
    float *       out4;
    double *      out8;
    long double * out16;
  };

  rho4 = (float *) field.values(id_field);
  out4 = (float*) initialize_output_(field.field_data());

  double vol = hx*hy*hz;
  
  switch (precision) {
  case precision_single:
    if (out4) {
#ifdef DEBUG_REFINE      
      double min =  std::numeric_limits<float>::max();
      double max = -min;
      double sum = 0.0;
#endif      
      for (int iz=gz; iz<mz-gz; iz++) {
	for (int iy=gy; iy<my-gy; iy++) {
	  for (int ix=gx; ix<mx-gx; ix++) {
	    int i = ix + mx*(iy + my*iz);
	    double mass = vol*rho4[i];
#ifdef DEBUG_REFINE      
	    min = std::min(min,mass);
	    max = std::max(max,mass);
	    sum += mass;
#endif	    
	    if      (mass < mass_max_coarsen) out4[i] = -1;
	    else if (mass < mass_min_refine)  out4[i] =  0;
	    else                              out4[i] = +1;
	  }
	}
      }
#ifdef DEBUG_REFINE      
      double sum2= 0.0;
      for (int i=0; i<mx*my*mz; i++) sum2 += vol*rho4[i];
      CkPrintf ("DEBUG_REFINE_MASS %s sum (%18.15g) [%18.15g]\n",
		name_.c_str(),sum,sum2);
      CkPrintf ("DEBUG_REFINE_MASS %s refine %18.15g\n",
		name_.c_str(),mass_min_refine);
#endif      
    }
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  double mass = vol*rho4[i];
	  if (mass > mass_min_refine)  any_refine  = true;
	  if (mass > mass_max_coarsen) all_coarsen = false;
	}
      }
    }
    break;
  case precision_double:
    if (out8) {
      double min =  std::numeric_limits<double>::max();
      double max = -min;
      for (int iz=gz; iz<mz-gz; iz++) {
	for (int iy=gy; iy<my-gy; iy++) {
	  for (int ix=gx; ix<mx-gx; ix++) {
	    int i = ix + mx*(iy + my*iz);
	    double mass = vol*rho8[i];
	    min = std::min(min,mass);
	    max = std::max(max,mass);
	    if      (mass < mass_max_coarsen) out8[i] = -1;
	    else if (mass < mass_min_refine)  out8[i] =  0;
	    else                              out8[i] = +1;
	  }
	}
      }
    }
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  double mass = vol*rho8[i];
	  if (mass > mass_min_refine)  any_refine  = true;
	  if (mass > mass_max_coarsen) all_coarsen = false;
	}
      }
    }
    break;
  case precision_quadruple:
    if (rho16) {
      for (int iz=gz; iz<mz-gz; iz++) {
	for (int iy=gy; iy<my-gy; iy++) {
	  for (int ix=gx; ix<mx-gx; ix++) {
	    int i = ix + mx*(iy + my*iz);
	    long double mass = vol*rho16[i];
	    if      (mass < mass_max_coarsen) rho16[i] = -1;
	    else if (mass < mass_min_refine)  rho16[i] =  0;
	    else                              rho16[i] = +1;
	  }
	}
      }
    }
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i = ix + mx*(iy + my*iz);
	  long double mass = vol*rho16[i];
	  if (mass > mass_min_refine)  any_refine  = true;
	  if (mass > mass_max_coarsen) all_coarsen = false;
	}
      }
    }
    break;
  default:
    ERROR2("EnzoRefineMass::apply",
	   "Unknown precision %d for field %d",
	   precision,0);
    break;
  }

  int adapt_result = 
    any_refine ?  adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

  // Don't refine if already at maximum level
  adjust_for_level_( &adapt_result, block->level() );

  return adapt_result;

}

//======================================================================

