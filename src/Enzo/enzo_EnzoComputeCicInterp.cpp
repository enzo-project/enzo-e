// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeCicInterp.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-05-05
/// @brief    Implements the EnzoComputeCicInterp class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputeCicInterp::EnzoComputeCicInterp  
(FieldDescr    * field_descr,
 std::string     field_name,
 ParticleDescr * particle_descr,
 std::string     particle_type,
 std::string     particle_attribute)
  : it_p_ (particle_descr->type_index (particle_type)),
    ia_p_ (particle_descr->attribute_index (it_p_,particle_attribute)),
    if_ (field_descr->field_id (field_name))
{
}

//----------------------------------------------------------------------

void EnzoComputeCicInterp::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);

  p | it_p_;
  p | ia_p_;
  p | if_;
}

//----------------------------------------------------------------------

void EnzoComputeCicInterp::compute ( Block * block) throw()
{

  if (!block->is_leaf()) return;

  Field    field    = block->data()->field();
  Particle particle = block->data()->particle();

  int field_precision    = field.precision(if_);
  int particle_precision = particle.attribute_type (it_p_,ia_p_);

  ASSERT1 ("EnzoComputeCicInterp::compute()",
	   "Unsupported particle_precision %d",
	   particle_precision,
	   particle_precision == type_single ||
	   particle_precision == type_double);

  ASSERT1 ("EnzoComputeCicInterp::compute()",
	   "Unsupported field_precision %d",
	   field_precision,
	   field_precision == precision_single ||
	   field_precision == precision_double);

  int p_single = (particle_precision == type_single);
  int p_double = (particle_precision == type_double);
  int f_single = (field_precision == precision_single);
  int f_double = (field_precision == precision_double);

  if (p_single && f_single) compute_<float, float>  (block); 
  if (p_single && f_double) compute_<float, double> (block); 
  if (p_double && f_single) compute_<double,float> (block); 
  if (p_double && f_double) compute_<double,double>(block); 

}

//----------------------------------------------------------------------

template <typename TP, typename TF>
void EnzoComputeCicInterp::compute_(Block * block)
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = enzo_block->data()->field();
  Particle particle = enzo_block->data()->particle();

  TF * vf = (TF*)field.values(if_);

  const int ia_x = particle.attribute_position(it_p_,0);
  const int ia_y = particle.attribute_position(it_p_,1);
  const int ia_z = particle.attribute_position(it_p_,2);

  const int dp =  particle.stride(it_p_,ia_x);
  const int da =  particle.stride(it_p_,ia_p_);

  const int rank = block->rank();

  int mx,my,mz;
  field.dimensions(0,&mx,&my,&mz);
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  // Get block extents and cell widths
  double xm,ym,zm;
  double xp,yp,zp;
  block->lower(&xm,&ym,&zm);
  block->upper(&xp,&yp,&zp);

  const int nb = particle.num_batches(it_p_);

  //  CkPrintf ("TRACE %s:%d rank = %d\n",__FILE__,__LINE__,rank);
  //  CkPrintf ("TRACE %s:%d nb = %d\n",__FILE__,__LINE__,nb);

  for (int ib=0; ib<nb; ib++) {

    TP * vp = (TP*) particle.attribute_array(it_p_, ia_p_, ib);

    const int np = particle.num_particles(it_p_,ib);
    //    CkPrintf ("TRACE %s:%d np = %d\n",__FILE__,__LINE__,np);

    if (rank == 1) {

      double * xa = (double *) particle.attribute_array (it_p_,ia_x,ib);

      for (int ip=0; ip<np; ip++) {

	double x = xa[ip*dp];

	double tx = nx*(x - xm) / (xp - xm) - 0.5;

	int ix0 = gx + floor(tx);

	int ix1 = ix0 + 1;

	double x0 = 1.0 - (tx - floor(tx));

	double x1 = 1.0 - x0;

	vp[ip*da] = x0*vf[ix0] + x1*vf[ix1];

      }

    } else if (rank == 2) {

      double * xa = (double *) particle.attribute_array (it_p_,ia_x,ib);
      double * ya = (double *) particle.attribute_array (it_p_,ia_y,ib);

      //      CkPrintf ("TRACE %s:%d\n",__FILE__,__LINE__);
      for (int ip=0; ip<np; ip++) {

	double x = xa[ip*dp];
	double y = ya[ip*dp];

	double tx = nx*(x - xm) / (xp - xm) - 0.5;
	double ty = ny*(y - ym) / (yp - ym) - 0.5;

	int ix0 = gx + floor(tx);
	int iy0 = gy + floor(ty);

	int ix1 = ix0 + 1;
	int iy1 = iy0 + 1;

	double x0 = 1.0 - (tx - floor(tx));
	double y0 = 1.0 - (ty - floor(ty));

	double x1 = 1.0 - x0;
	double y1 = 1.0 - y0;

	if ( ! ( (0.0 <= x0 && x0 <= 1.0) ||
		 (0.0 <= y0 && y0 <= 1.0) ||
		 (0.0 <= x1 && x1 <= 1.0) ||
		 (0.0 <= y1 && y1 <= 1.0))) {
	  CkPrintf ("ERROR? %s:%d [xy][01] = %f %f  %f %f\n",
		    __FILE__,__LINE__,x0,y0,x1,y1);
	}
	vp[ip*da] = x0*y0*vf[ix0+mx*iy0] 
	  +         x1*y0*vf[ix1+mx*iy0] 
	  +         x0*y1*vf[ix0+mx*iy1] 
	  +         x1*y1*vf[ix1+mx*iy1];

      }

    } else if (rank == 3) {

      double * xa = (double *) particle.attribute_array (it_p_,ia_x,ib);
      double * ya = (double *) particle.attribute_array (it_p_,ia_y,ib);
      double * za = (double *) particle.attribute_array (it_p_,ia_z,ib);

      for (int ip=0; ip<np; ip++) {

	double x = xa[ip*dp];
	double y = ya[ip*dp];
	double z = za[ip*dp];

	double tx = nx*(x - xm) / (xp - xm) - 0.5;
	double ty = ny*(y - ym) / (yp - ym) - 0.5;
	double tz = nz*(z - zm) / (zp - zm) - 0.5;

	int ix0 = gx + floor(tx);
	int iy0 = gy + floor(ty);
	int iz0 = gz + floor(tz);

	int ix1 = ix0 + 1;
	int iy1 = iy0 + 1;
	int iz1 = iz0 + 1;

	double x0 = 1.0 - (tx - floor(tx));
	double y0 = 1.0 - (ty - floor(ty));
	double z0 = 1.0 - (tz - floor(tz));

	double x1 = 1.0 - x0;
	double y1 = 1.0 - y0;
	double z1 = 1.0 - z0;

	vp[ip*da] = x0*y0*z0*vf[ix0+mx*(iy0+my*iz0)] 
	  +         x1*y0*z0*vf[ix1+mx*(iy0+my*iz0)] 
	  +         x0*y1*z0*vf[ix0+mx*(iy1+my*iz0)] 
	  +         x1*y1*z0*vf[ix1+mx*(iy1+my*iz0)]
	  +         x0*y0*z1*vf[ix0+mx*(iy0+my*iz1)] 
	  +         x1*y0*z1*vf[ix1+mx*(iy0+my*iz1)] 
	  +         x0*y1*z1*vf[ix0+mx*(iy1+my*iz1)] 
	  +         x1*y1*z1*vf[ix1+mx*(iy1+my*iz1)];

      }

    }
  }

}

