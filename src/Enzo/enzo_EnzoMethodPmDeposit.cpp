// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPmDeposit.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPmDeposit class
///
/// The EnzoMethodPmDeposit method computes a "density_total" field,
/// which includes the "density" field plus mass from gravitating
/// particles (particles in the "mass" group, e.g. "dark" matter
/// particles)

#include "cello.hpp"
#include "enzo.hpp"

// #define DEBUG_COLLAPSE

#define FORTRAN_NAME(NAME) NAME##_

extern "C" void  FORTRAN_NAME(dep_grid_cic)
  (enzo_float * de,enzo_float * de_t,enzo_float * temp,
   enzo_float * vx, enzo_float * vy, enzo_float * vz,
   enzo_float * dt, enzo_float * rfield, int *rank,
   enzo_float * hx, enzo_float * hy, enzo_float * hz,
   int * mx,int * my,int * mz,
   int * gxi,int * gyi,int * gzi,
   int * nxi,int * nyi,int * nzi,
   int * ,int * ,int * ,
   int * nx,int * ny,int * nz,
   int * ,int * ,int * );

//----------------------------------------------------------------------

EnzoMethodPmDeposit::EnzoMethodPmDeposit ( double alpha)
  : Method(),
    alpha_(alpha)
{
  const int rank = cello::rank();

  cello::define_field ("density");
  cello::define_field ("density_total");
  cello::define_field ("density_particle");
  cello::define_field ("density_particle_accumulate");
  if (rank >= 1) cello::define_field ("velocity_x");
  if (rank >= 2) cello::define_field ("velocity_y");
  if (rank >= 3) cello::define_field ("velocity_z");

  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());

  Refresh * refresh = cello::refresh(ir_post_);

  refresh->add_field("density");
  refresh->add_field("velocity_x");
  refresh->add_field("velocity_y");
  refresh->add_field("velocity_z");
}

//----------------------------------------------------------------------

void EnzoMethodPmDeposit::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | alpha_;
}

//----------------------------------------------------------------------

void EnzoMethodPmDeposit::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    Particle particle (block->data()->particle());
    Field    field    (block->data()->field());

    int rank = cello::rank();
    const int level = block->level();
    
    enzo_float * de_t = (enzo_float *)
      field.values("density_total");
    enzo_float * de_p = (enzo_float *)
      field.values("density_particle");
    enzo_float * de_pa = (enzo_float *)
      field.values("density_particle_accumulate");

    int mx,my,mz;
    field.dimensions(0,&mx,&my,&mz);
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);

    // NOTE: density_total is now cleared in EnzoMethodGravity to
    // instead of here to possible race conditions with refresh.  This
    // means EnzoMethodPmDeposit ("pm_deposit") currently CANNOT be
    // used without EnzoMethodGravity ("gravity")

    const int m = mx*my*mz;
    std::fill_n(de_p,m,0.0);
    std::fill_n(de_pa,m,0.0);
    
    // Get block extents and cell widths
    double xm,ym,zm;
    double xp,yp,zp;
    double hx,hy,hz;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);
    block->cell_width(&hx,&hy,&hz);

    // To calculate densities from particle masses, we need the inverse
    // volume of each cell
    double inv_vol = 1.0;
    if (rank >= 1) inv_vol /= hx;
    if (rank >= 2) inv_vol /= hy;
    if (rank >= 3) inv_vol /= hz;
    
    // Get cosmological scale factors, if cosmology is turned on
    enzo_float cosmo_a=1.0;
    enzo_float cosmo_dadt=0.0;
    EnzoPhysicsCosmology * cosmology = enzo::cosmology();
    
    if (cosmology) {
      cosmology->compute_expansion_factor(&cosmo_a,&cosmo_dadt,
					  block->time() + alpha_*block->dt());
    }

    const double dt = alpha_ * block->dt() / cosmo_a;

    // Get the number of particle types in the 'has_mass' group
    ParticleDescr * particle_descr = cello::particle_descr();
    Grouping * particle_groups = particle_descr->groups();
    const int num_mass = particle_groups->size("has_mass");

    enzo_float * pdens = NULL;
   
    // Loop over all particles that have mass
    // To do: rename has_mass to is_gravitating
    for (int ipt = 0; ipt < num_mass; ipt++) {
      const int it = particle.type_index(particle_groups->item("has_mass",ipt));

      // Index for mass / density attribute / constant
      int ind = 0;

      // If particle has constant called "mass" or "density", create an
      // array which will be later filled with a constant density value
      if (particle.has_constant(it, "mass") ||
	  particle.has_constant(it,"density"))
	pdens = new enzo_float[particle_descr->batch_size()];
      	
      // check correct precision for position
      int ia = particle.attribute_index(it,"x");
      int ba = particle.attribute_bytes(it,ia); // "bytes (actual)"
      int be = sizeof(enzo_float);                // "bytes (expected)"
      
      ASSERT4 ("EnzoMethodPmUpdate::compute()",
               "Particle type %s attribute %s defined as %s but expecting %s",
               particle.type_name(it).c_str(),
               particle.attribute_name(it,ia).c_str(),
               ((ba == 4) ? "single" : ((ba == 8) ? "double" : "quadruple")),
               ((be == 4) ? "single" : ((be == 8) ? "double" : "quadruple")),
               (ba == be));

 
      // Loop over batches
      for (int ib=0; ib<particle.num_batches(it); ib++) {
      
        const int np = particle.num_particles(it,ib);

	// First case: particle type has a constant called "mass"
	if (particle.has_constant (it,"mass")){

	  // Check if particle doesn't also have an attribute called "mass"
	  ASSERT1("EnzoMethodPmDeposit::compute",
		  "Particle type %s has both a constant and an attribute called"
		  "'mass'. Exiting.", particle.type_name(it).c_str(),
		  !particle.has_attribute(it,"mass"));

	  // Check if particle doesn't also have a constant called "density"
	  ASSERT1("EnzoMethodPmDeposit::compute",
		  "Particle type %s has both a constant called 'mass' and a "
		  "constant called 'density'. Exiting.",
		  particle.type_name(it).c_str(),
		  !particle.has_constant(it,"density"));

	  // Check if particle doesn't also have an attribute called "density"
	  ASSERT1("EnzoMethodPmDeposit::compute",
		  "Particle type %s has both a constant called 'mass' and an "
		  "attribute called 'density'. Exiting.",
		  particle.type_name(it).c_str(),
		  !particle.has_attribute(it,"density"));

	  // In this case we fill the first np elements of pdens with the constant
	  // mass value multiplied by inv_vol
	  ind = particle.constant_index(it,"mass");
	  for (int ip = 0; ip<np; ip++)
	    pdens[ip] = *((enzo_float *)(particle.constant_value (it,ind))) * inv_vol;
          
        } else if (particle.has_attribute(it,"mass")) {

	  // Check if particle doesn't also have a constant called "density"
	  ASSERT1("EnzoMethodPmDeposit::compute",
		  "Particle type %s has both an attribute called 'mass' and a "
		  "constant called 'density'. Exiting.",
		  particle.type_name(it).c_str(),
		  !particle.has_constant(it,"density"));

	  // Check if particle doesn't also have an attribute called "density"
	  ASSERT1("EnzoMethodPmDeposit::compute",
		  "Particle type %s has both an attribute called 'mass' and an "
		  "attribute called 'density'. Exiting.",
		  particle.type_name(it).c_str(),
		  !particle.has_attribute(it,"density"));

	  // In this case we set pdens to point to the mass attribute array
	  // and then multiply all its elements by inv_vol
	  ind = particle.attribute_index(it,"mass");
	  pdens = (enzo_float *) particle.attribute_array( it, ind, ib);
	  for (int ip = 0; ip<np; ip++)
	    pdens[ip] *= inv_vol;
        } else if (particle.has_constant(it,"density")) {

	  // Check if particle doesn't also have an attribute called "density"
	  ASSERT1("EnzoMethodPmDeposit::compute",
		  "Particle type %s has both a constant and an "
		  "attribute called 'density'. Exiting.",
		  particle.type_name(it).c_str(),
		  !particle.has_attribute(it,"density"));
	  
	  // In this case we fill the first np elements of pdens with the constant
	  // density value, and then rescale according to the refinement level
	  ind = particle.constant_index(it,"density");
	  for (int ip = 0; ip<np; ip++){
	    pdens[ip] = *((enzo_float *)(particle.constant_value (it,ind)));
	    for (int j = 0; j<level*rank; j++) pdens[ip] *= 2.0;
	  }
	} else if (particle.has_attribute(it,"density")) {
	  
	  // In this case we set pdens to point to the mass attribute array
	  // and then rescale the values according to the refinement level 
	  ind = particle.attribute_index(it,"density");
	  pdens = (enzo_float *) particle.attribute_array( it, ind, ib);
	  for (int ip = 0; ip<np ; ip++)
	    for (int j = 0 ; j<level*rank; j++)
	      pdens[ip] *= 2.0;
	  
	} else {

	  ERROR1("EnzoMethodPmDeposit::compute",
		"Particle type %s has neither a constant nor an attribute" 
                "called 'mass' or `density', yet is in the 'has_mass' group",
		 particle.type_name(it).c_str());
	}
	
	// If mass / density is a constant, we simply loop through the pdens
	// array, and so dm = 1. If its an attribute, we need to get
	// the stride length
	
	const int stride = (particle.has_attribute(it,"mass") ||
			    particle.has_attribute(it,"density"))
	                   ? particle.stride(it,ind) : 1;

	// Deposit densities to the grid with CIC scheme
        if (rank == 1) {

	  const int ia_x  = particle.attribute_index(it,"x");
	  const int ia_vx = particle.attribute_index(it,"vx");
	  
	  enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
	  enzo_float * vxa = (enzo_float *)particle.attribute_array (it,ia_vx,ib);

      	  const int dp =  particle.stride(it,ia_x);
      	  const int dv =  particle.stride(it,ia_vx);
	  
#ifdef DEBUG_COLLAPSE
          CkPrintf ("DEBUG_COLLAPSE vxa[0] = %lg\n",vxa[0]);
#endif

	  for (int ip=0; ip<np; ip++) {
	    
      	    double x = xa[ip*dp] + vxa[ip*dv]*dt;

      	    double tx = nx*(x - xm) / (xp - xm) - 0.5;
	    
	    int ix0 = gx + floor(tx);
	    
      	    int ix1 = ix0 + 1;
	    
      	    double x0 = 1.0 - (tx - floor(tx));
      	    double x1 = 1.0 - x0;
	    
            de_p[ix0] += pdens[ip*stride]*x0;
            de_p[ix1] += pdens[ip*stride]*x1;

	    if (de_p[ix0] < 0.0) {
      	      CkPrintf ("%s:%d ERROR: de_p %d = %f\n",
      	   	        __FILE__,__LINE__,ix0,de_p[ix0]);
      	    }

	    if (de_p[ix1] < 0.0) {
      	      CkPrintf ("%s:%d ERROR: de_p %d = %f\n",
      	   	        __FILE__,__LINE__,ix0,de_p[ix1]);
      	    }
	    
      	  } // np
	  
        } else if (rank == 2) {

      	  const int ia_x  = particle.attribute_index(it,"x");
      	  const int ia_y  = particle.attribute_index(it,"y");
      	  const int ia_vx = particle.attribute_index(it,"vx");
      	  const int ia_vy = particle.attribute_index(it,"vy");
	  
      	  // Batch arrays
      	  enzo_float * xa  = (enzo_float *)particle.attribute_array (it,ia_x,ib);
      	  enzo_float * ya  = (enzo_float *)particle.attribute_array (it,ia_y,ib);
      	  enzo_float * vxa = (enzo_float *)particle.attribute_array (it,ia_vx,ib);
      	  enzo_float * vya = (enzo_float *)particle.attribute_array (it,ia_vy,ib);

      	  const int dp =  particle.stride(it,ia_x);
      	  const int dv =  particle.stride(it,ia_vx);
	  
	  for (int ip=0; ip<np; ip++) {
	    
      	    double x = xa[ip*dp] + vxa[ip*dv]*dt;
      	    double y = ya[ip*dp] + vya[ip*dv]*dt;

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
	    
      	    de_p[ix0+mx*iy0] += pdens[ip*stride]*x0*y0;
      	    de_p[ix1+mx*iy0] += pdens[ip*stride]*x1*y0;
      	    de_p[ix0+mx*iy1] += pdens[ip*stride]*x0*y1;
      	    de_p[ix1+mx*iy1] += pdens[ip*stride]*x1*y1;
	    
      	    if (de_p[ix0+mx*iy0] < 0.0) {
      	      CkPrintf ("%s:%d ERROR: de_p %d %d = %f\n",
      	   	        __FILE__,__LINE__,ix0,iy0,de_p[ix0+mx*iy0]);
      	    }
      	    if (de_p[ix1+mx*iy0] < 0.0) {
      	      CkPrintf ("%s:%d ERROR: de_p %d %d = %f\n",
      		        __FILE__,__LINE__,ix1,iy0,de_p[ix1+mx*iy0]);
      	    }
      	    if (de_p[ix0+mx*iy1] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d = %f\n",
      		        __FILE__,__LINE__,ix0,iy1,de_p[ix0+mx*iy1]);
      	    }
      	    if (de_p[ix1+mx*iy1] < 0.0) {
      	      CkPrintf ("%s:%d ERROR: de_p %d %d = %f\n",
      		        __FILE__,__LINE__,ix1,iy1,de_p[ix1+mx*iy1]);
      	    }
	  } // ip

        } else if (rank == 3) {

      	  const int ia_x  = particle.attribute_index(it,"x");
      	  const int ia_y  = particle.attribute_index(it,"y");
      	  const int ia_z  = particle.attribute_index(it,"z");
      	  const int ia_vx = particle.attribute_index(it,"vx");
      	  const int ia_vy = particle.attribute_index(it,"vy");
      	  const int ia_vz = particle.attribute_index(it,"vz");
	  
      	  enzo_float * xa  = (enzo_float *) particle.attribute_array (it,ia_x,ib);
      	  enzo_float * ya  = (enzo_float *) particle.attribute_array (it,ia_y,ib);
      	  enzo_float * za  = (enzo_float *) particle.attribute_array (it,ia_z,ib);

      	  // Particle batch velocities
      	  enzo_float * vxa = (enzo_float *) particle.attribute_array (it,ia_vx,ib);
      	  enzo_float * vya = (enzo_float *) particle.attribute_array (it,ia_vy,ib);
      	  enzo_float * vza = (enzo_float *) particle.attribute_array (it,ia_vz,ib);

#ifdef DEBUG_COLLAPSE
          CkPrintf ("DEBUG_COLLAPSE vxa[0] = %lg\n",vxa[0]);
#endif

      	  const int dp =  particle.stride(it,ia_x);
      	  const int dv =  particle.stride(it,ia_vx);

	  for (int ip=0; ip<np; ip++) {
	    
	    // Copy batch particle velocities to temporary block field velocities
	    
	    double x = xa[ip*dp] + vxa[ip*dv]*dt;
	    double y = ya[ip*dp] + vya[ip*dv]*dt;
	    double z = za[ip*dp] + vza[ip*dv]*dt;
	    
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
	    
            de_p[ix0+mx*(iy0+my*iz0)] += pdens[ip*stride]*x0*y0*z0;
            de_p[ix1+mx*(iy0+my*iz0)] += pdens[ip*stride]*x1*y0*z0;
            de_p[ix0+mx*(iy1+my*iz0)] += pdens[ip*stride]*x0*y1*z0;
            de_p[ix1+mx*(iy1+my*iz0)] += pdens[ip*stride]*x1*y1*z0;
            de_p[ix0+mx*(iy0+my*iz1)] += pdens[ip*stride]*x0*y0*z1;
            de_p[ix1+mx*(iy0+my*iz1)] += pdens[ip*stride]*x1*y0*z1;
            de_p[ix0+mx*(iy1+my*iz1)] += pdens[ip*stride]*x0*y1*z1;
            de_p[ix1+mx*(iy1+my*iz1)] += pdens[ip*stride]*x1*y1*z1;

	    if (de_p[ix0+mx*(iy0+my*iz0)] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d %d = %f\n",
			__FILE__,__LINE__,
			ix0,iy0,iz0,de_p[ix0+mx*(iy0+my*iz0)]);
	    }
	    if (de_p[ix1+mx*(iy0+my*iz0)] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d %d = %f\n",
			__FILE__,__LINE__,
			ix1,iy0,iz0,de_p[ix1+mx*(iy0+my*iz0)]);
	    }
	    if (de_p[ix0+mx*(iy1+my*iz0)] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d %d = %f\n",
			__FILE__,__LINE__,
			ix0,iy1,iz0,de_p[ix0+mx*(iy1+my*iz0)]);
	    }
	    if (de_p[ix1+mx*(iy1+my*iz0)] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d %d = %f\n",
			__FILE__,__LINE__,
			ix1,iy1,iz0,de_p[ix1+mx*(iy1+my*iz0)]);
	    }
	    if (de_p[ix0+mx*(iy0+my*iz1)] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d %d = %f\n",
			__FILE__,__LINE__,
			ix0,iy0,iz1,de_p[ix0+mx*(iy0+my*iz1)]);
	    }
	    if (de_p[ix1+mx*(iy0+my*iz1)] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d %d = %f\n",
			__FILE__,__LINE__,
			ix1,iy0,iz1,de_p[ix1+mx*(iy0+my*iz1)]);
	    }
	    if (de_p[ix0+mx*(iy1+my*iz1)] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d %d = %f\n",
			__FILE__,__LINE__,
			ix0,iy1,iz1,de_p[ix0+mx*(iy1+my*iz1)]);
	    }
	    if (de_p[ix1+mx*(iy1+my*iz1)] < 0.0) {
	      CkPrintf ("%s:%d ERROR: de_p %d %d %d = %f\n",
			__FILE__,__LINE__,
			ix1,iy1,iz1,de_p[ix1+mx*(iy1+my*iz1)]);
	    }
          } // Loop over particle in batch
        } // if rank == 3
	
      } // Loop over batches
      
      // If constant mass, delete the pdens array
      if (particle.has_constant(it,"mass") ||
	  particle.has_constant(it,"density"))
	delete [] pdens;
      
    } // Loop over particle types in 'has_mass' group
    
    //--------------------------------------------------
    // Add gas density
    //--------------------------------------------------
    enzo_float * de = (enzo_float *) field.values("density");

    enzo_float * temp =   new enzo_float [4*m];
    enzo_float * de_gas = new enzo_float [m];
    enzo_float * rfield = new enzo_float [m];

    std::fill_n(temp, 4*m, 0.0);
    std::fill_n(de_gas, m, 0.0);
    std::fill_n(rfield, m, 0.0);

    int gxi=gx;
    int gyi=gy;
    int gzi=gz;
    int nxi=mx-gx-1;
    int nyi=my-gy-1;
    int nzi=mz-gz-1;
    int i0 = 0;
    int i1 = 1;

    /* Stefan: Not sure if these next four lines are correct.
       I include the factors of cosmo_a for consistency with 
       previous version. Also, not sure why dtf is
       initialised with the value of alpha_.
     */
    
    enzo_float hxf = hx * cosmo_a;
    enzo_float hyf = hy * cosmo_a;
    enzo_float hzf = hz * cosmo_a;
    enzo_float dtf = alpha_;

    enzo_float * vxf = (enzo_float *) field.values("velocity_x");
    enzo_float * vyf = (enzo_float *) field.values("velocity_y");
    enzo_float * vzf = (enzo_float *) field.values("velocity_z");

    enzo_float * vx = new enzo_float [m];
    enzo_float * vy = new enzo_float [m];
    enzo_float * vz = new enzo_float [m];

    for (int i=0; i<m; i++) vx[i] = vxf[i];

    if (rank >= 2) for (int i=0; i<m; i++) vy[i] = vyf[i];
    else           for (int i=0; i<m; i++) vy[i] = 0.0;

    if (rank >= 3) for (int i=0; i<m; i++) vz[i] = vzf[i];
    else           for (int i=0; i<m; i++) vz[i] = 0.0;

    
    FORTRAN_NAME(dep_grid_cic)(de,de_gas,temp,
			       vx, vy, vz,
			       &dtf, rfield, &rank,
			       &hxf,&hyf,&hzf,
			       &mx,&my,&mz,
			       &gxi,&gyi,&gzi,
			       &nxi,&nyi,&nzi,
			       &i0,&i0,&i0,
			       &nx,&ny,&nz,
			       &i1,&i1,&i1);

    for (int i=0; i<mx*my*mz; i++) {
      de_t[i] = de_pa[i] = de_p[i];
    }

    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
    	  for (int ix=gx; ix<mx-gx; ix++) {
	        int i = ix + mx*(iy + my*iz);
	        int ig = (ix-gx) + nx*((iy-gy) + ny*(iz-gz));
	        de_t[i] += de_gas[ig];
    	  }
      }
    }

    delete [] rfield;
    delete [] temp;
    delete [] vx;
    delete [] vy;
    delete [] vz;

    delete [] de_gas;
  }

  block->compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodPmDeposit::timestep ( Block * block ) throw()
{
  double dt = std::numeric_limits<double>::max();

  return dt;
}
