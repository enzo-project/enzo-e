// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPmUpdate.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPmUpdate class

#include "cello.hpp"
#include "enzo.hpp"

// #define DEBUG_PM_UPDATE

#ifdef DEBUG_PM_UPDATE
#  define TRACE_PM(MESSAGE)						\
  CkPrintf ("%s:%d %s\n",						\
	    __FILE__,__LINE__,MESSAGE);				
#else
#  define TRACE_PM(MESSAGE) /* ... */
#endif

//----------------------------------------------------------------------

EnzoMethodPmUpdate::EnzoMethodPmUpdate 
( const FieldDescr * field_descr,
  const ParticleDescr * particle_descr, double max_dt ) 
  : Method(),
    max_dt_(max_dt)
{
  TRACE_PM("EnzoMethodPmUpdate()");
  // Initialize default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_particles(particle_descr->num_types());
  refresh(ir)->add_all_fields(field_descr->field_count());

  // PM parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodPmUpdate::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | max_dt_;
}

//----------------------------------------------------------------------

void EnzoMethodPmUpdate::compute ( Block * block) throw()
{
  TRACE_PM("compute()");

  if (block->is_leaf()) {

    const int rank = block->rank();

    FieldDescr    * fd = block->data()->field_descr();
    ParticleDescr * pd = block->data()->particle_descr();

    if (rank >= 1) {
      EnzoComputeCicInterp interp_x (fd, "acceleration_x", pd, "dark", "ax");
      interp_x.compute(block);
    }

    if (rank >= 2) {
      EnzoComputeCicInterp interp_y (fd, "acceleration_y", pd, "dark", "ay");
      interp_y.compute(block);
    }

    if (rank >= 3) {
      EnzoComputeCicInterp interp_z (fd, "acceleration_z", pd, "dark", "az");
      interp_z.compute(block);
    }

    Particle particle = block->data()->particle();

    const int it = particle.type_index ("dark");

    const int ia_x  = (rank >= 1) ? particle.attribute_index (it, "x") : -1;
    const int ia_y  = (rank >= 2) ? particle.attribute_index (it, "y") : -1;
    const int ia_z  = (rank >= 3) ? particle.attribute_index (it, "z") : -1;

    const int ia_vx = (rank >= 1) ? particle.attribute_index (it, "vx") : -1;
    const int ia_vy = (rank >= 2) ? particle.attribute_index (it, "vy") : -1;
    const int ia_vz = (rank >= 3) ? particle.attribute_index (it, "vz") : -1;

    const int ia_ax = (rank >= 1) ? particle.attribute_index (it, "ax") : -1;
    const int ia_ay = (rank >= 2) ? particle.attribute_index (it, "ay") : -1;
    const int ia_az = (rank >= 3) ? particle.attribute_index (it, "az") : -1;

    const int dp = particle.stride(it, ia_x);
    const int dv = particle.stride(it, ia_vx);
    const int da = particle.stride(it, ia_ax);

    const int nb = particle.num_batches (it);

    const double dt = block->dt();

    for (int ib=0; ib<nb; ib++) {

      double *x=0, *y=0, *z=0;
      double *vx=0, *vy=0, *vz=0;
      double *ax=0, *ay=0, *az=0;

      if (rank >= 1) {
	x  = (double *) particle.attribute_array (it, ia_x,  ib);
	vx = (double *) particle.attribute_array (it, ia_vx, ib);
	ax = (double *) particle.attribute_array (it, ia_ax, ib);
      }
      if (rank >= 2) {
	y  = (double *) particle.attribute_array (it, ia_y,  ib);
	vy = (double *) particle.attribute_array (it, ia_vy, ib);
	ay = (double *) particle.attribute_array (it, ia_ay, ib);
      }
      if (rank >= 3) {
	z  = (double *) particle.attribute_array (it, ia_z,  ib);
	vz = (double *) particle.attribute_array (it, ia_vz, ib);
	az = (double *) particle.attribute_array (it, ia_az, ib);
      }

      const int np = particle.num_particles(it,ib);

      if (rank >= 1) {

	for (int ip=0; ip<np; ip++) {

	  const int ipdv = ip*dv;
	  const int ipdp = ip*dp;
	  const int ipda = ip*da;

	  vx[ipdv] += ax[ipda]*dt/2;
	  x [ipdp] += vx[ipdv]*dt;
	  vx[ipdv] += ax[ipda]*dt/2;

	}
	
      }
      if (rank >= 2) {

	for (int ip=0; ip<np; ip++) {

	  const int ipdv = ip*dv;
	  const int ipdp = ip*dp;
	  const int ipda = ip*da;

	  vy[ipdv] += ay[ipda]*dt/2;
	  y [ipdp] += vy[ipdv]*dt;
	  vy[ipdv] += ay[ipda]*dt/2;

	}
	
      }
      if (rank >= 3) {

	for (int ip=0; ip<np; ip++) {

	  const int ipdv = ip*dv;
	  const int ipdp = ip*dp;
	  const int ipda = ip*da;

	  vz[ipdv] += az[ipda]*dt/2;
	  z [ipdp] += vz[ipdv]*dt;
	  vz[ipdv] += az[ipda]*dt/2;

	}
      }
    }
  }

  block->compute_done(); 
  
}

//----------------------------------------------------------------------

double EnzoMethodPmUpdate::timestep ( Block * block ) const throw()
{
  TRACE_PM("timestep()");

  const int rank = block->rank();

  double dt = std::numeric_limits<double>::max();

  if (block->is_leaf()) {

    Particle particle = block->data()->particle();
    Field    field    = block->data()->field();

    const int it = particle.type_index ("dark");
    const int nb = particle.num_batches (it);

    const int ia_vx = particle.attribute_index (it, "vx");
    const int ia_vy = particle.attribute_index (it, "vy");
    const int ia_vz = particle.attribute_index (it, "vz");

    const int ia_ax = particle.attribute_index (it, "ax");
    const int ia_ay = particle.attribute_index (it, "ay");
    const int ia_az = particle.attribute_index (it, "az");

    const int dv = particle.stride(it, ia_vx);
    const int da = particle.stride(it, ia_ax);

    double xm,ym,zm;
    double xp,yp,zp;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);

    int nx,ny,nz;
    field.size(&nx,&ny,&nz);

    // ENZO particle time step
    //
    // for (dim = 0; dim < GridRank; dim++) {
    //   float dCell = CellWidth[dim][0]*a;
    //   for (i = 0; i < NumberOfParticles; i++) {
    //     dtTemp = dCell/max(fabs(ParticleVelocity[dim][i]), tiny_number);
    // 	dtParticles = min(dtParticles, dtTemp);
    //   }
    // }
 
    // /* Multiply resulting dt by ParticleCourantSafetyNumber. */
 
    // dtParticles *= ParticleCourantSafetyNumber;

    const double hx = (xp-xm)/nx;
    const double hy = (yp-ym)/ny;
    const double hz = (zp-zm)/nz;
    
    for (int ib=0; ib<nb; ib++) {
      const int np = particle.num_particles(it,ib);

      if (rank >= 1) {
	const double * vx = (const double *) 
	  particle.attribute_array (it, ia_vx, ib); 
	const double * ax = (const double *) 
	  particle.attribute_array (it, ia_ax, ib); 
	for (int ip=0; ip<np; ip++) {
	  const double v = fabs(vx[ip*dv]);
	  const double a = fabs(ax[ip*da]);
	  const double dt_v = hx /MAX(v,1e-6);
	  const double dt_a = sqrt(2.0*hx/MAX(a,1e-6));
	  dt = MIN(dt,dt_v);
	  dt = MIN(dt,dt_a);
	}
      }

      if (rank >= 2) {
	const double * vy = (const double *) 
	  particle.attribute_array (it, ia_vy, ib); 
	const double * ay = (const double *) 
	  particle.attribute_array (it, ia_ay, ib); 
	for (int ip=0; ip<np; ip++) {
	  const double v = fabs(vy[ip*dv]);
	  const double a = fabs(ay[ip*da]);
	  const double dt_v = hy /MAX(v,1e-6);
	  const double dt_a = sqrt(2.0*hy/MAX(a,1e-6));
	  dt = MIN(dt,dt_v);
	  dt = MIN(dt,dt_a);
	}
      }

      if (rank >= 3) {
	const double * vz = (const double *) 
	  particle.attribute_array (it, ia_vz, ib); 
	const double * az = (const double *) 
	  particle.attribute_array (it, ia_az, ib); 
	for (int ip=0; ip<np; ip++) {
	  const double v = fabs(vz[ip*dv]);
	  const double a = fabs(az[ip*da]);
	  const double dt_v = hz/MAX(v,1e-6);
	  const double dt_a = sqrt(2.0*hz/MAX(a,1e-6));
	  dt = MIN(dt,dt_v);
	  dt = MIN(dt,dt_a);
	}
      }

    }
  }
  dt = MIN(dt,max_dt_);
  return dt;
}
