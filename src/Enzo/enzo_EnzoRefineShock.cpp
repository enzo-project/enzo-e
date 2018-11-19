// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRefineShock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 21 16:02:39 PDT 2014
/// @brief    Implementation of EnzoRefineShock class

#include "enzo.hpp"

// #define DEBUG_ENZO_REFINE_SHOCK

//----------------------------------------------------------------------

EnzoRefineShock::EnzoRefineShock(double pressure_min_refine,
				 double pressure_max_coarsen,
				 double energy_ratio_min_refine,
				 double energy_ratio_max_coarsen,
				 double gamma,
				 bool comoving_coordinates,
				 int max_level,
				 bool include_ghosts,
				 std::string output) throw ()
  : Refine (0.0, 0.0, max_level, include_ghosts, output),
    pressure_min_refine_ (pressure_min_refine),
    pressure_max_coarsen_(pressure_max_coarsen),
    energy_ratio_min_refine_ (energy_ratio_min_refine),
    energy_ratio_max_coarsen_(energy_ratio_max_coarsen),
    gamma_(gamma),
    comoving_coordinates_(comoving_coordinates)
{
}

//----------------------------------------------------------------------

int EnzoRefineShock::apply ( Block * block ) throw ()
{

  Field field = block->data()->field();

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  int rank = cello::rank();

  // compute pressure using the EnzoComputePressure class

  EnzoComputePressure compute_pressure (gamma_,comoving_coordinates_);
  compute_pressure.compute(block);

  Data * data = block->data();
  double xm[3],xp[3];
  data->lower(&xm[0],&xm[1],&xm[2]);
  data->upper(&xp[0],&xp[1],&xp[2]);

  int id_velocity = field.field_id("velocity_x");

  ASSERT("EnzoRefineShock::apply",
	  "velocity_x field must be defined",
	 (id_velocity >= 0));

    void * v3[3] = {
    (rank >= 1) ? field.values("velocity_x") : NULL,
    (rank >= 2) ? field.values("velocity_y") : NULL,
    (rank >= 3) ? field.values("velocity_z") : NULL
  };

  void * te = field.values("total_energy");
  void * de = field.values("density");
  void * p  = field.values("pressure");
   
  int gx,gy,gz;
  field.ghost_depth(id_velocity, &gx,&gy,&gz);

  int nxd = nx + 2*gx;
  int nyd = ny + 2*gy;
  int nzd = nz + 2*gz;

  if (rank < 1) gx = 0;
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  //  precision_type precision = field.precision(id_velocity);

  void * output = initialize_output_(field.field_data());

  bool any_refine  = false;
  bool all_coarsen = true;
  
  evaluate_block_((const enzo_float**) v3,
		  (const enzo_float*)  te,
		  (const enzo_float*)  de,
		  (const enzo_float*)  p,
		  (enzo_float*) output,
		  nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		  &any_refine,&all_coarsen, rank);

  int adapt_result =  
    any_refine ? adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

  // Don't refine if already at maximum level
  adjust_for_level_ (&adapt_result,block->level());

  return adapt_result;

}

//----------------------------------------------------------------------

void EnzoRefineShock::evaluate_block_
(const enzo_float * v3[],
 const enzo_float * te,
 const enzo_float * de,
 const enzo_float * p,
 enzo_float * output,
 int ndx, int ndy, int ndz,
 int nx, int ny, int nz,
 int gx, int gy, int gz,
 bool *any_refine,
 bool * all_coarsen, 
 int rank)
{
  (*all_coarsen) = true;
  (*any_refine)  = false;

  const int d3[3] = {1, ndx, ndx*ndy};

#ifdef DEBUG_ENZO_REFINE_SHOCK
  enzo_float dp_min = std::numeric_limits<enzo_float>::max();
  enzo_float dp_max = -std::numeric_limits<enzo_float>::max();
  enzo_float er_min = std::numeric_limits<enzo_float>::max();
  enzo_float er_max = -std::numeric_limits<enzo_float>::max();
#endif
  
  for (int axis=0; axis<rank; axis++) {

    for (int iz=gz; iz<nz+gz; iz++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int ix=gx; ix<nx+gx; ix++) {

	  int i = ix + ndx*(iy + ndy*iz);
	  int id = d3[axis];

	  enzo_float dp = fabs    (p[i+id] - p[i-id]) 
	    / (std::min(p[i+id] , p[i-id])) ;

	  enzo_float dv = v3[axis][i+id] - v3[axis][i-id];

	  enzo_float e = p[i]/(gamma_ - 1.0);

	  enzo_float ep = te[i+id]*de[i+id];
	  enzo_float e0 = te[i]   *de[i];
	  enzo_float em = te[i-id]*de[i-id];

	  enzo_float er = e / std::max (std::max(em,e0),ep);

	  bool l_refine = (dv < 0.0) && 
	    (dp > pressure_min_refine_) &&
	    (er > energy_ratio_min_refine_);

	  bool l_same = (dv < 0.0) &&
	    (dp > pressure_max_coarsen_) &&
	    (er > energy_ratio_max_coarsen_);

#ifdef DEBUG_ENZO_REFINE_SHOCK
	  dp_min = std::min(dp_min,dp);
	  dp_max = std::max(dp_max,dp);
	  er_min = std::min(er_min,er);
	  er_max = std::max(er_max,er);
#endif
	  if (l_refine)  *any_refine = true;
	  if (l_same)    *all_coarsen = false;

	  if (output) {
	    if (l_same)   output[i] =  0;
	    if (l_refine) output[i] = +1;
	  }
	}
      }
    }
  }
#ifdef DEBUG_ENZO_REFINE_SHOCK
  CkPrintf ("%s dp limits %lf %lf  dp min/max %lf %lf\n",
	    pressure_max_coarsen_,pressure_min_refine_,
	    dp_min,dp_max);
  CkPrintf ("%s er limits %lf %lf  er min/max %lf %lf\n",
	    energy_ratio_max_coarsen_,energy_ratio_min_refine_,
	    er_min,er_max);
#endif  
}
//======================================================================

