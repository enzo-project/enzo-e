// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodBackgroundAcceleration.cpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     2018-05
/// @brief    Implements the EnzoMethodBackgroundAcceleration class


#include "cello.hpp"
#include "enzo.hpp"

// What are these defs and do I need them?
//#include "enzo.decl.h"
//
//
// #define CK_TEMPLATES_ONLY
// #include "enzo.def.h"
// #undef CK_TEMPLATES_ONLY

//---------------------------------------------------------------------

EnzoMethodBackgroundAcceleration::EnzoMethodBackgroundAcceleration
(bool zero_acceleration) // do need
 : Method(),
   zero_acceleration_(zero_acceleration),
   mx_(0), my_(0), mz_(0),
   gx_(0), gy_(0), gz_(0),
   xm_(0), ym_(0), zm_(0),
   hx_(0), hy_(0), hz_(0)
{

  this->G_four_pi_ = 4.0 * cello::pi * cello::grav_constant;

  FieldDescr * field_descr = cello::field_descr();

  //const int id =  field_descr->field_id("density");
  const int iax = field_descr->field_id("acceleration_x");
  const int iay = field_descr->field_id("acceleration_y");
  const int iaz = field_descr->field_id("acceleration_z");


  // Do not need to refresh acceleration fields in this method
  // since we do not need to know any ghost zone information
  //
  const int ir = add_refresh(4,0,neighbor_leaf,sync_neighbor,
                             enzo_sync_id_method_background_acceleration);

  refresh(ir)->add_field(iax);
  refresh(ir)->add_field(iay);
  refresh(ir)->add_field(iaz);

  // refresh(ir)->add_field(id);
  // iax iay and iax are used in method gravity to do
  // refreshing of fields.... Do I need ot do this here? If so
  // why do I need it..... otherwise I can just worry about field
  // access with field_values?

  return;

} // EnzoMethodBackgroundAcceleration

//-------------------------------------------------------------------
void EnzoMethodBackgroundAcceleration::compute ( Block * block) throw()
{
  if (block->is_leaf()){
    this->compute_(block);
  }

  // Wait for all blocks to finish before continuing
  block->compute_done();
  return;
}

void EnzoMethodBackgroundAcceleration::compute_ (Block * block) throw()
{
  /* This applies background potential to grid and particles */

  //TRACE_METHOD("compute()",block);
  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  if (!(enzo_config->method_background_acceleration_apply_acceleration)) return;

  Field field = block->data()->field();

  // Obtain grid sizes and ghost sizes

  field.dimensions (0,&mx_,&my_,&mz_);
  field.ghost_depth(0,&gx_,&gy_,&gz_);
  double xp, yp, zp;
  block->data()->lower(&xm_,&ym_,&zm_);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm_,xp,&hx_,ym_,yp,&hy_,zm_,zp,&hz_);

  // We will probably never be in the situation of constant acceleration
  // and cosmology, but just in case.....
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  enzo_float cosmo_a = 1.0;

  const int rank = cello::rank();

  if (cosmology) {
    enzo_float cosmo_dadt = 0.0;
    double dt    = block->dt();
    double time  = block->time();
    cosmology->compute_expansion_factor(&cosmo_a,&cosmo_dadt,time+0.5*dt);
    if (rank >= 1) hx_ *= cosmo_a;
    if (rank >= 2) hy_ *= cosmo_a;
    if (rank >= 3) hz_ *= cosmo_a;
  }


  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");

  int m = mx_ * my_ * mz_;
  if (zero_acceleration_){
    if (ax){ for(int i = 0; i < m; i ++){ ax[i] = 0.0;}}
    if (ay){ for(int i = 0; i < m; i ++){ ay[i] = 0.0;}}
    if (az){ for(int i = 0; i < m; i ++){ az[i] = 0.0;}}
  }

  Particle particle = enzo_block->data()->particle();

  if (enzo_config->method_background_acceleration_type == "GalaxyModel"){

    this->GalaxyModel(ax, ay, az, &particle, rank,
                      cosmo_a, enzo_config, enzo_units, enzo_block->dt);

  } else if (enzo_config->method_background_acceleration_type == "PointMass"){
    this->PointMass(ax, ay, az, &particle, rank,
                    cosmo_a, enzo_config, enzo_units, enzo_block->dt);
  } else {

    ERROR("EnzoMethodBackgroundAcceleration::compute_()",
          "Background acceleration type not recognized");

  }


  return;

} // compute

void EnzoMethodBackgroundAcceleration::PointMass(enzo_float * ax,
                                                 enzo_float * ay,
                                                 enzo_float * az,
                                                 Particle * particle,
                                                 const int rank,
                                                 const enzo_float cosmo_a,
                                                 const EnzoConfig * enzo_config,
                                                 const EnzoUnits * enzo_units,
                                                 const double dt)
                                                 throw() {

  // just need to define position of each cell

  double mass = enzo_config->method_background_acceleration_mass *
                cello::mass_solar / enzo_units->mass();
  double rcore = std::max(0.1*hx_,
                     enzo_config->method_background_acceleration_core_radius/enzo_units->length());
  double G = this->G_four_pi_ *
            enzo_units->density() * enzo_units->time() * enzo_units->time();

  double x = 0.0, y = 0.0, z = 0.0;

  for (int iz=0; iz<mz_; iz++){
    if (rank >= 3) z = zm_ + (iz - gz_ + 0.5)*hz_ - enzo_config->method_background_acceleration_center[2];

    for (int iy=0; iy<my_; iy++){
      if (rank >= 2) y = ym_ + (iy - gy_ + 0.5)*hy_ - enzo_config->method_background_acceleration_center[1];

      for (int ix=0; ix<mx_; ix++){
        x = xm_ + (ix - gx_ + 0.5)*hx_ - enzo_config->method_background_acceleration_center[0];

        double rsqr  = x*x + y*y + z*z;
        double r     = sqrt(rsqr);

        double accel = G * std::min(mass / ((rsqr)*r*cosmo_a),
                             mass / ((rcore*rcore*rcore)*cosmo_a));

        int i = INDEX(ix,iy,iz,mx_,my_);

        if (ax) ax[i] -= accel * x;
        if (ay) ay[i] -= accel * y;
        if (az) az[i] -= accel * z;

      }
    }
  } // end loop over grid cells

  // Update particle accelerations here -- leave for now

  // Particle particle = block->data()->particle();

  return;
}

void EnzoMethodBackgroundAcceleration::GalaxyModel(enzo_float * ax,
                                                   enzo_float * ay,
                                                   enzo_float * az,
                                                   Particle * particle,
                                                   const int rank,
                                                   const enzo_float cosmo_a,
                                                   const EnzoConfig * enzo_config,
                                                   const EnzoUnits * enzo_units,
                                                   const double dt)
                                                   throw() {

  double DM_mass     = enzo_config->method_background_acceleration_DM_mass *
                         cello::mass_solar / enzo_units->mass();
  double DM_mass_radius = enzo_config->method_background_acceleration_DM_mass_radius *
                          cello::kpc_cm / enzo_units->length();
  double DM_density  = enzo_config->method_background_acceleration_DM_density /
                         enzo_units->density();
  double stellar_r   = enzo_config->method_background_acceleration_stellar_scale_height_r *
                         cello::kpc_cm / enzo_units->length();
  double stellar_z   = enzo_config->method_background_acceleration_stellar_scale_height_z *
                         cello::kpc_cm / enzo_units->length();
  double stellar_mass = enzo_config->method_background_acceleration_stellar_mass *
                        cello::mass_solar / enzo_units->mass();
  double bulge_mass   = enzo_config->method_background_acceleration_bulge_mass *
                        cello::mass_solar / enzo_units->mass();
  double bulgeradius = enzo_config->method_background_acceleration_bulge_radius *
                        cello::kpc_cm / enzo_units->length();
  const double * amom = enzo_config->method_background_acceleration_angular_momentum;

  double G = this->G_four_pi_ *
             enzo_units->density() * enzo_units->time() * enzo_units->time();
  double G_code = cello::grav_constant * enzo_units->density() * enzo_units->time() * enzo_units->time();

  double rcore = enzo_config->method_background_acceleration_core_radius *
                 cello::kpc_cm / enzo_units->length();

  if (DM_mass > 0.0){
    double xtemp = DM_mass_radius / rcore;

    // compute the density constant for an NFW halo (rho_o)
    DM_density = (DM_mass / (4.0 * cello::pi * std::pow(rcore,3))) /
                       (std::log(1.0+xtemp)-xtemp/(1.0+xtemp));

  } else {
    double xtemp = DM_mass_radius / rcore;

    DM_mass = 4.0 * cello::pi / 3.0 * (std::pow(rcore,3) * DM_density) *
                 3.0 * (std::log(1.0 + xtemp) - xtemp/(1.0+xtemp));
  }


  double x = 0.0, y = 0.0, z = 0.0;

  double accel_sph, accel_R, accel_z;

  for (int iz=0; iz<mz_; iz++){
     if (rank >= 3) z = zm_ + (iz - gz_ + 0.5)*hz_ - enzo_config->method_background_acceleration_center[2];

     for (int iy=0; iy<my_; iy++){
       if (rank >= 2) y = ym_ + (iy - gy_ + 0.5)*hy_ - enzo_config->method_background_acceleration_center[1];

       for (int ix=0; ix<mx_; ix++){
         x = xm_ + (ix - gx_ + 0.5)*hx_ - enzo_config->method_background_acceleration_center[0];

         // double rsqr  = x*x + y*y + z*z;
         // double r     = sqrt(rsqr);

         double zheight = amom[0]*x + amom[1]*y + amom[2]*z; // height above disk

         // projected positions in plane of the disk
         double xplane = x - zheight*amom[0];
         double yplane = y - zheight*amom[1];
         double zplane = z - zheight*amom[2];

         double radius = sqrt(xplane*xplane + yplane*yplane + zplane*zplane + zheight*zheight);
         double rcyl   = sqrt(xplane*xplane + yplane*yplane + zplane*zplane);

         // need to multiple all of the below by the gravitational constants
         double xtemp     = radius/rcore;
         double Rtemp     = DM_mass_radius / rcore;

         //double
         accel_sph = G * bulge_mass / pow(radius + bulgeradius,2) +    // bulge
                     + 4.0 * G_code * cello::pi * DM_density * rcore *
                          (log(1.0+xtemp) - (xtemp / (1.0+xtemp))) / (xtemp*xtemp);

         //double
         accel_R   = G * stellar_mass * rcyl / sqrt( pow( pow(rcyl,2)
                             + pow(stellar_r + sqrt( pow(zheight,2)
                            + pow(stellar_z,2)),2),3));
         //double
         accel_z   = G * stellar_mass / sqrt(pow(zheight,2)
                               + pow(stellar_z,2))*zheight/sqrt(pow(pow(rcyl,2)
                               + pow(stellar_r + sqrt(pow(zheight,2)
                               + pow(stellar_z,2)),2),3))
                                   * (stellar_z * sqrt(pow(zheight,2) + pow(stellar_z,2)));

         accel_sph = (radius  == 0.0 ? 0.0 : std::fabs(accel_sph) / (radius*cosmo_a));
         accel_R   = (rcyl    == 0.0 ? 0.0 : std::fabs(accel_R)   / (rcyl*cosmo_a));
         accel_z   = (zheight == 0.0 ? 0.0 : std::fabs(accel_z)*zheight/std::fabs(zheight) / cosmo_a);

         accel_R = 0.0; accel_z = 0.0;
         // now apply accelerations in cartesian (grid) coordinates
         int i = INDEX(ix,iy,iz,mx_,my_);
         if (ax) ax[i] -= (accel_sph * x + accel_R*xplane + accel_z*amom[0]);
         if (ay) ay[i] -= (accel_sph * y + accel_R*yplane + accel_z*amom[1]);
         if (az) az[i] -= (accel_sph * z + accel_R*zplane + accel_z*amom[2]);
        }
     }
  } // end loop over grid cells

  // Update particle accelerations for particles with mass

  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();

  int num_mass = particle_groups->size("has_mass");

  double dt_shift = 0.5 * dt;

  // Loop through particles to apply this to
  for (int ipt = 0; ipt < num_mass; ipt++){

    std::string particle_type = particle_groups->item("has_mass",ipt);
    int it = particle->type_index(particle_type);

    if (particle->num_particles(it) > 0){

      int ia_x = (rank >= 1) ? particle->attribute_index (it, "x") : -1;
      int ia_y = (rank >= 2) ? particle->attribute_index (it, "y") : -1;
      int ia_z = (rank >= 3) ? particle->attribute_index (it, "z") : -1;

      int ia_ax = (rank >= 1) ? particle->attribute_index (it, "ax") : -1;
      int ia_ay = (rank >= 2) ? particle->attribute_index (it, "ay") : -1;
      int ia_az = (rank >= 3) ? particle->attribute_index (it, "az") : -1;

      int dp = particle->stride(it, ia_x);
      int da = particle->stride(it, ia_ax);

      int nb = particle->num_batches (it);

      for (int ib=0; ib<nb; ib++){
        enzo_float *px=0, *py=0, *pz=0;
        enzo_float *pax=0, *pay=0, *paz=0;

        px  = (enzo_float *) particle->attribute_array (it, ia_x, ib);
        pax = (enzo_float *) particle->attribute_array (it, ia_ax, ib);
        py  = (enzo_float *) particle->attribute_array (it, ia_y, ib);
        pay = (enzo_float *) particle->attribute_array (it, ia_ay, ib);
        pz  = (enzo_float *) particle->attribute_array (it, ia_z, ib);
        paz = (enzo_float *) particle->attribute_array (it, ia_az, ib);

        int np = particle->num_particles(it,ib);

        for (int ip = 0; ip<np; ip++){
          int ipdp = ip*dp;
          int ipda = ip*da;

          x = px[ipdp] - enzo_config->method_background_acceleration_center[0];
          y = py[ipdp] - enzo_config->method_background_acceleration_center[1];
          z = pz[ipdp] - enzo_config->method_background_acceleration_center[2];

          double zheight = amom[0]*x + amom[1]*y + amom[2]*z; // height above disk

          // projected positions in plane of the disk
          double xplane = x - zheight*amom[0];
          double yplane = y - zheight*amom[1];
          double zplane = z - zheight*amom[2];

          double radius = sqrt(xplane*xplane + yplane*yplane + zplane*zplane + zheight*zheight);
          double rcyl   = sqrt(xplane*xplane + yplane*yplane + zplane*zplane);

          // need to multiple all of the below by the gravitational constants
          double xtemp     = radius/rcore;
          double Rtemp     = DM_mass_radius / rcore;

          //double
          accel_sph = G * bulge_mass / pow(radius + bulgeradius,2) +    // bulge
                     + 4.0 * G_code * cello::pi * DM_density * rcore *
                          (log(1.0+xtemp) - (xtemp / (1.0+xtemp))) / (xtemp*xtemp);

          accel_R   = G * stellar_mass * rcyl / sqrt( pow( pow(rcyl,2)
                              + pow(stellar_r + sqrt( pow(zheight,2)
                             + pow(stellar_z,2)),2),3));
          //double
          accel_z   = G * stellar_mass / sqrt(pow(zheight,2)
                                + pow(stellar_z,2))*zheight/sqrt(pow(pow(rcyl,2)
                                + pow(stellar_r + sqrt(pow(zheight,2)
                                + pow(stellar_z,2)),2),3))
                                    * (stellar_z * sqrt(pow(zheight,2) + pow(stellar_z,2)));

          accel_sph = (radius  == 0.0 ? 0.0 : std::fabs(accel_sph) / (radius*cosmo_a));
          accel_R   = (rcyl    == 0.0 ? 0.0 : std::fabs(accel_R)   / (rcyl*cosmo_a));
          accel_z   = (zheight == 0.0 ? 0.0 : std::fabs(accel_z)*zheight/std::fabs(zheight) / cosmo_a);

          // now apply accelerations in cartesian (grid) coordinates
          if (pax) pax[ipda] -= (accel_sph * x + accel_R*xplane + accel_z*amom[0]);
          if (pay) pay[ipda] -= (accel_sph * y + accel_R*yplane + accel_z*amom[1]);
          if (paz) paz[ipda] -= (accel_sph * z + accel_R*zplane + accel_z*amom[2]);


        } // end loop over particles

      } // end loop over batch

    } // end if particles exist

  } // end loop over particles


  return;
}

double EnzoMethodBackgroundAcceleration::timestep (Block * block) const throw()
{
  // Use the same timestep check as implemented for gravity. This
  // just goes through the acceleration fields and checkes to make
  // sure they
  Field field = block->data()->field();
  int mx, my, mz;
  int gx, gy, gz;
  field.dimensions(0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);


  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");

  enzo_float dt = std::numeric_limits<enzo_float>::max();

  double hx,hy,hz;
  block->cell_width(&hx,&hy,&hz);

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  const int rank = cello::rank();
  if (cosmology) {
    enzo_float cosmo_a = 1.0;
    enzo_float cosmo_dadt = 0.0;
    double dt   = block->dt();
    double time = block->time();
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);
    if (rank >= 1) hx*=cosmo_a;
    if (rank >= 2) hy*=cosmo_a;
    if (rank >= 3) hz*=cosmo_a;
  }

  if (ax) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hx/(fabs(ax[i]+1e-20)))));
	}
      }
    }
  }
  if (ay) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hy/(fabs(ay[i]+1e-20)))));
	}
      }
    }
  }
  if (az) {
    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
	for (int ix=gx; ix<mx-gx; ix++) {
	  int i=ix + mx*(iy + iz*my);
	  dt = std::min(enzo_float(dt),enzo_float(sqrt(hz/(fabs(az[i]+1e-20)))));
	}
      }
    }
  }

  return 0.5*dt;
}

/* Placeholder for the general DM potential
   methods (NFW is a specific case)
void EnzoMethodBackgroundAcceleration::General(void){

  return;
}
*/
