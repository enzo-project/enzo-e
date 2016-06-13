// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineParticleCount.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-05-19
/// @brief    Implementation of RefineParticleCount class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineParticleCount::RefineParticleCount
(
 double min_refine,
 double max_coarsen,
 int    max_level,
 bool include_ghosts,
 std::string output,
 double level_exponent
) throw ()
  : Refine(min_refine,
	   max_coarsen,
	   max_level,
	   include_ghosts,
	   output),
    level_exponent_(level_exponent)
{
  TRACE("RefineParticleCount::RefineParticleCount");
}

//----------------------------------------------------------------------

int RefineParticleCount::apply ( Block * block ) throw ()
{

  Particle particle = block->data()->particle();
  Field    field    = block->data()->field();

  const int level = block->level();

  const int rank = block->rank();
  double xm,ym,zm;
  double xp,yp,zp;
  block->lower (&xm,&ym,&zm);
  block->upper (&xp,&yp,&zp);

  const int it = particle.type_index ("dark");

  const int ia_x = particle.attribute_position (it,0);
  const int ia_y = particle.attribute_position (it,1);
  const int ia_z = particle.attribute_position (it,2);

  const int nb = particle.num_batches (it);
  const int dp = particle.stride(it,ia_x);

  size_t count = 0;

  for (int ib=0; ib<nb; ib++) {

    const int np = particle.num_particles (it,ib);

    if (rank == 1) {

      double * xa = (double *) particle.attribute_array (it,ia_x,ib);

      for (int ip=0; ip<np; ip++) {
	const double x = xa[ip*dp];
	if (xm <= x && x <= xp) ++count;
      }

    } else if (rank == 2) {

      double * xa = (double *) particle.attribute_array (it,ia_x,ib);
      double * ya = (double *) particle.attribute_array (it,ia_y,ib);

      for (int ip=0; ip<np; ip++) {
	double x = xa[ip*dp];
	double y = ya[ip*dp];
	if ((xm <= x && x <= xp) &&
	    (ym <= y && y <= yp)) ++count;
      }

    } else if (rank == 3) {

      double * xa = (double *) particle.attribute_array (it,ia_x,ib);
      double * ya = (double *) particle.attribute_array (it,ia_y,ib);
      double * za = (double *) particle.attribute_array (it,ia_z,ib);

      for (int ip=0; ip<np; ip++) {
	double x = xa[ip*dp];
	double y = ya[ip*dp];
	double z = za[ip*dp];

	if ((xm <= x && x <= xp) &&
	    (ym <= y && y <= yp) &&
	    (zm <= z && z <= zp)) ++count;

      }
    }
  }

  int adapt_result = 
    (count > min_refine_)  ? adapt_refine :
    ((count < max_coarsen_) ? adapt_coarsen : adapt_same);

  // Don't refine if already at maximum level
  adjust_for_level_( &adapt_result, level );

  return adapt_result;

}


//======================================================================

