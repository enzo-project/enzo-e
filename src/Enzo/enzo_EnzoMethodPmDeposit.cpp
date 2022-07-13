// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPmDeposit.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Stefan Arridge (stefan.arridge@gmail.com)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPmDeposit class
///
/// The EnzoMethodPmDeposit method computes a "density_total" field,
/// which includes the "density" field plus mass from gravitating
/// particles (particles in the "is_gravitating" group, e.g. "dark" matter
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
  // Check if particle types in "is_gravitating" group have either a constant
  // or an attribute called "mass" (but not both).
  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();
  const int num_is_grav = particle_groups->size("is_gravitating");
  for (int ipt = 0; ipt < num_is_grav; ipt++) {
    const int it = particle_descr->type_index(particle_groups->item("is_gravitating",ipt));
    
    // Count number of attributes or constants called "mass",
    // which should be equal to 1
    int num_mass = 0;
    if (particle_descr->has_constant (it,"mass")) ++num_mass;
    if (particle_descr->has_attribute (it,"mass")) ++num_mass;

    ASSERT1("EnzoMethodPmDeposit::EnzoMethodPmDeposit",
	    "Particle type %s, in the \"is_gravitating\" group, "
            "must have either an attribute or a constant "
	    "called \"mass\" (but not both) . Exiting.",
	    particle_descr->type_name(it).c_str(),
	    num_mass == 1);
  }
  
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

namespace { // define local helper functions in anonymous namespace

  /// deposits mass from all gravitating particles onto density_particle_arr
  ///
  /// @param[out] density_particle_arr The array where the deposited mass
  ///     density is stored
  /// @param[in]  Block Contains the particle data to use for accumulation
  /// @param[in]  dt_div_cosmoa Length of time to "drift" the particles before
  ///     before deposition divided by the scale factor (computed for the time
  ///     after particles have been drifted)
  /// @param[in]  inv_vol One divided by the volume of a cell. When
  ///     `cello::rank()` is 2, this is really "inverse-area" and when it's 1,
  ///     this is really "inverse-cell-width". In cosmological simulations this
  ///     should be the comoving quantity.
  /// @param[in]      mx,my,mz Specifies the number of cells along each
  ///     dimension of an array (including ghost cells)
  /// @param[in]      gx,gy,gz Specifies the number of cells in the ghost zone
  ///     for each dimensions
  void deposit_particles_(const CelloArray<enzo_float,3>& density_particle_arr,
                          Block* block, double dt_div_cosmoa, double inv_vol,
                          int mx, int my, int mz,
                          int gx, int gy, int gz)
  {
    Particle particle (block->data()->particle());
    Field    field    (block->data()->field());

    int rank = cello::rank();

    // compute extent of the active zone
    const int nx = mx - 2 * gx;
    const int ny = (rank >=2) ? my - 2 * gy : 1;
    const int nz = (rank >=3) ? mz - 2 * gz : 1;

    enzo_float * de_p = density_particle_arr.data();

    // Get block extents
    double xm,ym,zm;
    double xp,yp,zp;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);

    // Get the number of particle types in the "is_gravitating" group
    ParticleDescr * particle_descr = cello::particle_descr();
    Grouping * particle_groups = particle_descr->groups();
    const int num_is_grav = particle_groups->size("is_gravitating");

    // For particle types where "mass" is an attribute,
    // pmass will be set to point to array of particle masses
    // For particle types where "mass" is a constant,
    // pmass will point to the constant value.
    enzo_float * pmass = NULL;

    // The for the mass "array" (if "mass" is a constant, then
    // there won't be a mass array, and the stride will be set
    // to zero.
    int dm;

    // Loop over particle types in "is_gravitating" group
    for (int ipt = 0; ipt < num_is_grav; ipt++) {
      const int it = particle.type_index(particle_groups->item("is_gravitating",ipt));

      // Index for mass attribute / constant
      int imass = 0;

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

	if (particle.has_attribute(it,"mass")) {

	  // Particle type has an attribute called "mass".
	  // In this case we set pmass to point to the mass attribute array
	  // Also set dm to be the stride for the "mass" attribute
	  imass = particle.attribute_index(it,"mass");
	  pmass = (enzo_float *) particle.attribute_array( it, imass, ib);
	  dm = particle.stride(it,imass);

	} else {

	  // Particle type has a constant called "mass".
	  // In this case we set pmass to point to the value
	  // of the mass constant.
	  // dm is set to 0, which will mean that we can loop through an
	  // "array" of length 1.
	  imass = particle.constant_index(it,"mass");
	  pmass = (enzo_float*)particle.constant_value(it,imass);
	  dm = 0;
	}

	// Deposit densities to the grid with CIC scheme
	if (rank == 1) {

	  const int ia_x  = particle.attribute_index(it,"x");
	  const int ia_vx = particle.attribute_index(it,"vx");

	  enzo_float * xa =  (enzo_float *)particle.attribute_array (it,ia_x,ib);
	  enzo_float * vxa = (enzo_float *)particle.attribute_array (it,ia_vx,ib);
	  const int dp =  particle.stride(it,ia_x);
	  const int dv =  particle.stride(it,ia_vx);
#ifdef DEBUG_COLLAPS
	  CkPrintf ("DEBUG_COLLAPSE vxa[0] = %lg\n",vxa[0]);
#endif

	  for (int ip=0; ip<np; ip++) {
	    double x = xa[ip*dp] + vxa[ip*dv]*dt_div_cosmoa;

	    double tx = nx*(x - xm) / (xp - xm) - 0.5;

	    int ix0 = gx + floor(tx);
	    int ix1 = ix0 + 1;
	    double x0 = 1.0 - (tx - floor(tx));
	    double x1 = 1.0 - x0;

	    // Density is mass times inverse volume
	    // If mass is a constant, then dm is 0, pmass[ip * dm] is pmass[0], which
	    // just dereferences pmass.
	    enzo_float pdens = pmass[ip*dm] * inv_vol;
	    de_p[ix0] += pdens * x0;
	    de_p[ix1] += pdens * x1;

	    if (de_p[ix0] < 0.0)
	      WARNING3("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d] = %g",
		       block->name().c_str(),ix0,de_p[ix0]);

	    if (de_p[ix1] < 0.0)
	      WARNING3("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d] = %g",
		       block->name().c_str(),ix1,de_p[ix1]);

	  } // Loop over particles in batch

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

	    double x = xa[ip*dp] + vxa[ip*dv]*dt_div_cosmoa;
	    double y = ya[ip*dp] + vya[ip*dv]*dt_div_cosmoa;

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

	    // Density is mass times inverse volume
	    // If mass is a constant, then dm is 0, pmass[ip * dm] is pmass[0], which
	    // just dereferences pmass.
	    enzo_float pdens = pmass[ip*dm] * inv_vol;
	    de_p[ix0+mx*iy0] += pdens * x0 * y0;
	    de_p[ix1+mx*iy0] += pdens * x1 * y0;
	    de_p[ix0+mx*iy1] += pdens * x0 * y1;
	    de_p[ix1+mx*iy1] += pdens * x1 * y1;

	    if (de_p[ix0+mx*iy0] < 0.0)
	      WARNING4("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d] = %g",
		       block->name().c_str(),ix0,iy0,de_p[ix0+mx*iy0]);

	    if (de_p[ix1+mx*iy0] < 0.0)
	      WARNING4("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d] = %g",
		       block->name().c_str(),ix1,iy0,de_p[ix1+mx*iy0]);

	    if (de_p[ix0+mx*iy1] < 0.0)
	      WARNING4("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d] = %g",
		       block->name().c_str(),ix0,iy1,de_p[ix0+mx*iy1]);

	    if (de_p[ix1+mx*iy1] < 0.0)
	      WARNING4("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d] = %g",
		       block->name().c_str(),ix1,iy1,de_p[ix1+mx*iy1]);


	  } // Loop over particles in batch

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

	    double x = xa[ip*dp] + vxa[ip*dv]*dt_div_cosmoa;
	    double y = ya[ip*dp] + vya[ip*dv]*dt_div_cosmoa;
	    double z = za[ip*dp] + vza[ip*dv]*dt_div_cosmoa;

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

	    // Density is mass times inverse volume
	    // If mass is a constant, then dm is 0, pmass[ip * dm] is pmass[0], which
	    // just dereferences pmass.
	    enzo_float pdens = pmass[ip*dm] * inv_vol;
	    de_p[ix0+mx*(iy0+my*iz0)] += pdens * x0 * y0 * z0;
	    de_p[ix1+mx*(iy0+my*iz0)] += pdens * x1 * y0 * z0;
	    de_p[ix0+mx*(iy1+my*iz0)] += pdens * x0 * y1 * z0;
	    de_p[ix1+mx*(iy1+my*iz0)] += pdens * x1 * y1 * z0;
	    de_p[ix0+mx*(iy0+my*iz1)] += pdens * x0 * y0 * z1;
	    de_p[ix1+mx*(iy0+my*iz1)] += pdens * x1 * y0 * z1;
	    de_p[ix0+mx*(iy1+my*iz1)] += pdens * x0 * y1 * z1;
	    de_p[ix1+mx*(iy1+my*iz1)] += pdens * x1 * y1 * z1;

	    if (de_p[ix0+mx*(iy0+my*iz0)] < 0.0)
	      WARNING5("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d,%d] = %g",
		       block->name().c_str(),ix0,iy0,iz0,
		       de_p[ix0+mx*(iy0+my*iz0)]);

	    if (de_p[ix1+mx*(iy0+my*iz0)] < 0.0)
	      WARNING5("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d,%d] = %g",
		       block->name().c_str(),ix1,iy0,iz0,
		       de_p[ix1+mx*(iy0+my*iz0)]);

	    if (de_p[ix0+mx*(iy1+my*iz0)] < 0.0)
	      WARNING5("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d,%d] = %g",
		       block->name().c_str(),ix0,iy1,iz0,
		       de_p[ix0+mx*(iy1+my*iz0)]);

	    if (de_p[ix1+mx*(iy1+my*iz0)] < 0.0)
	      WARNING5("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d,%d] = %g",
		       block->name().c_str(),ix1,iy1,iz0,
		       de_p[ix1+mx*(iy1+my*iz0)]);

	    if (de_p[ix0+mx*(iy0+my*iz1)] < 0.0)
	      WARNING5("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d,%d] = %g",
		       block->name().c_str(),ix0,iy0,iz1,
		       de_p[ix0+mx*(iy0+my*iz1)]);

	    if (de_p[ix1+mx*(iy0+my*iz1)] < 0.0)
	      WARNING5("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d,%d] = %g",
		       block->name().c_str(),ix1,iy0,iz1,
		       de_p[ix1+mx*(iy0+my*iz1)]);

	    if (de_p[ix0+mx*(iy1+my*iz1)] < 0.0)
	      WARNING5("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d,%d] = %g",
		       block->name().c_str(),ix0,iy1,iz1,
		       de_p[ix0+mx*(iy1+my*iz1)]);

	    if (de_p[ix1+mx*(iy1+my*iz1)] < 0.0)
	      WARNING5("EnzoMethodPmDeposit",
		       "Block %s: de_p[%d,%d,%d] = %g",
		       block->name().c_str(),ix1,iy1,iz1,
		       de_p[ix1+mx*(iy1+my*iz1)]);

	  } // Loop over particles in batch
	} // if rank == 3

      } // Loop over batches

    } // Loop over particle types in "is_gravitating" group

  }

  //----------------------------------------------------------------------

  /// deposits mas density from gas onto density_tot_arr
  ///
  /// @param[in, out] density_tot_arr The array where density gets accumulated
  /// @param[in]      field Contains the field data to use for accumulation
  /// @param[in]      dt Length of time to "drift" the density field before
  ///     deposition divided by the scale factor at the deposition time
  /// @param[in]      hx_prop,hy_prop,hz_prop The width of cell along
  ///     each axis. These specify the proper lengths at the time that we
  ///     deposit the density (after any drift).
  /// @param[in]      mx,my,mz Specifies the number of cells along each
  ///     dimension of an array (including ghost cells)
  /// @param[in]      gx,gy,gz Specifies the number of cells in the ghost zone
  ///     for each dimensions
  void deposit_gas_(const CelloArray<enzo_float, 3>& density_tot_arr,
                    Field& field, enzo_float dt_div_cosmoa,
                    enzo_float hx_prop, enzo_float hy_prop, enzo_float hz_prop,
                    int mx, int my, int mz,
                    int gx, int gy, int gz){

    int rank = cello::rank();
    const int m = mx*my*mz;

    // compute extent of the active zone
    int nx = mx - 2 * gx;
    int ny = (rank >=2) ? my - 2 * gy : 1;
    int nz = (rank >=3) ? mz - 2 * gz : 1;

    // retrieve primary fields needed for depositing gas density
    enzo_float * de = (enzo_float *) field.values("density");
    enzo_float * vxf = (enzo_float *) field.values("velocity_x");
    enzo_float * vyf = (enzo_float *) field.values("velocity_y");
    enzo_float * vzf = (enzo_float *) field.values("velocity_z");

    // allocate and zero-initialize scratch arrays for missing velocity
    // components.
    std::vector<enzo_float> vel_scratch(m*(3 - rank), 0.0);
    if (rank < 2) vyf = vel_scratch.data();
    if (rank < 3) vzf = vel_scratch.data() + m;

    // allocate memory for storing the deposited gas density. The treatment of
    // this array is a little odd:
    //    - After calling dep_grid_cic, we treat this memory as though just
    //      enough space has been allocated for holding data in the active
    //      zone (as though it holds just nz*ny*nx entries). Note: this is
    //      mainly known from earlier implemenations of this code.
    //    - It seems as though this array MUST have an allocation for extra
    //      array allocations. If we allocate just nz*ny*nx entries, then
    //      dep_grid_cic ends up writing to invalid memory locations. Earlier
    //      versions of the code allocate mz*my*mx entries which seems
    //      sufficient.
    std::vector<enzo_float> deposited_gas_density(m, 0.0);

    // allocate temporary arrays
    std::vector<enzo_float> temp(4*m, 0.0); // scratch space
    std::vector<enzo_float> rfield(m, 0.0); // legacy argument relevant for
                                            // patch-based AMR

    // specify the indices corresponding to the start and end (inclusive) of
    // the source fields (namely de, vxf, vyf, vzf)
    int src_starti_x=gx;
    int src_starti_y=gy;
    int src_starti_z=gz;
    int src_endi_x=mx-gx-1;
    int src_endi_y=my-gy-1;
    int src_endi_z=mz-gz-1;

    // dummy args
    int offset_val = 0;
    int refine_factor = 1;

    FORTRAN_NAME(dep_grid_cic)(de, deposited_gas_density.data(), temp.data(),
                               vxf, vyf, vzf,
                               &dt_div_cosmoa, rfield.data(), &rank,
                               &hx_prop,&hy_prop,&hz_prop,
                               &mx,&my,&mz,
                               &src_starti_x, &src_starti_y, &src_starti_z,
                               &src_endi_x, &src_endi_y, &src_endi_z,
                               &offset_val, &offset_val, &offset_val,
                               &nx,&ny,&nz,
                               &refine_factor, &refine_factor, &refine_factor);

    // construct a subarray density_tot that just includes the active zone
    CelloArray<enzo_float,3> density_tot_active_zone = density_tot_arr.subarray
      (CSlice(gz, mz - gz), CSlice(gy, my - gy), CSlice(gx, mx - gx));

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
        for (int ix=0; ix<nx; ix++) {
          // treat deposited_gas_density like it's a contiguous 3D array with a
          // shape of (nz, ny, nx) that extends over the active zone
          int i_gas = ix + nx*(iy + ny*iz);
          density_tot_active_zone(iz,iy,ix) += deposited_gas_density[i_gas];
        }
      }
    }
  }

}

//----------------------------------------------------------------------

void EnzoMethodPmDeposit::compute ( Block * block) throw()
{

   if (enzo::simulation()->cycle() == enzo::config()->initial_cycle) {
    // Check if the gravity method is being used and that pm_deposit
    // precedes the gravity method.
    ASSERT("EnzoMethodPmDeposit",
           "Error: pm_deposit method must precede gravity method.",
           enzo::problem()->method_precedes("pm_deposit", "gravity"));
  }

  if (block->is_leaf()) {

    Field    field    (block->data()->field());

    CelloArray<enzo_float,3> density_tot_arr =
      field.view<enzo_float>("density_total");
    CelloArray<enzo_float,3> density_particle_arr =
      field.view<enzo_float>("density_particle");
    CelloArray<enzo_float,3> density_particle_accum_arr =
      field.view<enzo_float>("density_particle_accumulate");

    int mx,my,mz;
    field.dimensions(0,&mx,&my,&mz);
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);

    std::fill_n(density_particle_arr.data(), mx*my*mz, 0.0);

    // NOTE 2022-06-24: previously, we filled density_particle_accum_arr with
    // zeros at around this location and included the following note:
    //     NOTE: density_total is now cleared in EnzoMethodGravity to
    //     instead of here to possible race conditions with refresh.  This
    //     means EnzoMethodPmDeposit ("pm_deposit") currently CANNOT be
    //     used without EnzoMethodGravity ("gravity")
    // This operation & comment didn't sense since we completely overwrite
    // values of density_total & density_particle_accum_arr later in this method

    // Get cosmological scale factors, if cosmology is turned on
    enzo_float cosmo_a=1.0;
    enzo_float cosmo_dadt=0.0;
    EnzoPhysicsCosmology * cosmology = enzo::cosmology();
    if (cosmology) {
      cosmology->compute_expansion_factor(&cosmo_a,&cosmo_dadt,
                                          block->time() + alpha_*block->dt());
    }
    const double dt_div_cosmoa = alpha_ * block->dt() / cosmo_a;

    double hx,hy,hz;
    block->cell_width(&hx,&hy,&hz);

    // add mass from particles
    {
      int rank = cello::rank();
      // To calculate densities from particles with "mass" attributes
      // or constants, we need the inverse volume of cells in this block.
      double inv_vol = 1.0;
      if (rank >= 1) inv_vol /= hx;
      if (rank >= 2) inv_vol /= hy;
      if (rank >= 3) inv_vol /= hz;

      deposit_particles_(density_particle_arr, block, dt_div_cosmoa, inv_vol,
                         mx, my, mz,
                         gx, gy, gz);

      // update density_tot_arr
      density_particle_arr.copy_to(density_tot_arr);
      density_particle_arr.copy_to(density_particle_accum_arr);
    }


    // add mass from gas
    {
      // Grid_DepositBaryons.C from enzo-dev overrides the drift time for the
      // gas density to be zero when using the PPM and Zeus solvers.
      //
      // The way that the Gravity source terms are currently included with
      // EnzoMethodPpm and EnzoMethodMHDVlct is currently very similar (in
      // terms of temporal integration), so we also set this to zero.
      const double gas_dt_div_cosmoa = 0.0;

      // The use of "proper" cell-widths was carried over for consistency with
      // earlier versions of the code. Based on Grid_DepositBaryons.C from
      // enzo-dev, it seems like this may not be correct.
      deposit_gas_(density_tot_arr, field, gas_dt_div_cosmoa,
                   hx*cosmo_a, hy*cosmo_a, hz*cosmo_a,
                   mx, my, mz,
                   gx, gy, gz);
    }
  }

  block->compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodPmDeposit::timestep ( Block * block ) throw()
{
  double dt = std::numeric_limits<double>::max();

  return dt;
}
