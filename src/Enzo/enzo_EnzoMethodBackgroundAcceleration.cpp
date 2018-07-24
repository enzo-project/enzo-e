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

extern CProxy_EnzoSimulation proxy_enzo_simulation;

//---------------------------------------------------------------------

EnzoMethodBackgroundAcceleration::EnzoMethodBackgroundAcceleration
(const FieldDescr * field_descr, bool zero_acceleration) // do need
 : Method(),
   zero_acceleration_(zero_acceleration),
   mx_(0), my_(0), mz_(0),
   gx_(0), gy_(0), gz_(0),
   xm_(0), ym_(0), zm_(0),
   hx_(0), hy_(0), hz_(0)
{

  this->G_four_pi_ = 4.0 * cello::pi * cello::grav_constant;

  //const int id =  field_descr->field_id("density");
  //const int iax = field_descr->field_id("acceleration_x");
  //const int iay = field_descr->field_id("acceleration_y");
  //const int iaz = field_descr->field_id("acceleration_z");

/*
  // Do not need to refresh acceleration fields in this method
  // since we do not need to know any ghost zone information
  //
  const int ir = add_refresh(4,0,neighbor_leaf,sync_neighbor,
                             enzo_sync_id_method_background_acceleration);

  refresh(ir)->add_field(iax);
  refresh(ir)->add_field(iay);
  refresh(ir)->add_field(iaz);
*/
  //refresh(ir)->add_field(id);
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
  /* This applies background potential to grid. AE: need to write similar function
     for particles. Unsure if this should be a compute(particle) kind of deal
     or can just be tacked into here */

  //TRACE_METHOD("compute()",block);
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  const EnzoConfig * enzo_config = static_cast<const EnzoConfig*>
       (enzo_block->simulation()->config());
  EnzoSimulation * simulation = proxy_enzo_simulation.ckLocalBranch();
  EnzoUnits * units = (EnzoUnits *) simulation->problem()->units();
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
  EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology * )
    block->simulation()->problem()->physics("cosmology");
  enzo_float cosmo_a = 1.0;
  if (cosmology) {
    const int rank = block->rank();
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

  if (enzo_config->method_background_acceleration_type == "GalaxyModel"){

    this->GalaxyModel(ax, ay, az, block->rank(),
                      cosmo_a, enzo_config, units);

  } else if (enzo_config->method_background_acceleration_type == "PointMass"){
    this->PointMass(ax, ay, az, block->rank(),
                    cosmo_a, enzo_config, units);
  }

  return;

} // compute

void EnzoMethodBackgroundAcceleration::PointMass(enzo_float * ax,
                                                 enzo_float * ay,
                                                 enzo_float * az,
                                                 const int rank,
                                                 const enzo_float cosmo_a,
                                                 const EnzoConfig * enzo_config,
                                                 const EnzoUnits * units)
                                                 throw() {

  // just need to define position of each cell

  double mass = enzo_config->method_background_acceleration_mass *
                cello::mass_solar / units->mass();
  double rcore = std::max(0.1*hx_,
                     enzo_config->method_background_acceleration_core_radius/units->length());
  double G = this->G_four_pi_ *
            units->density() * units->time() * units->time();

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
                                                   const int rank,
                                                   const enzo_float cosmo_a,
                                                   const EnzoConfig * enzo_config,
                                                   const EnzoUnits * units)
                                                   throw() {

  double DM_mass     = enzo_config->method_background_acceleration_DM_mass *
                         cello::mass_solar / units->mass();
  double DM_mass_radius = enzo_config->method_background_acceleration_DM_mass_radius *
                          cello::kpc / units->length();
  double DM_density  = enzo_config->method_background_acceleration_DM_density /
                         units->density();
  double stellar_r   = enzo_config->method_background_acceleration_stellar_scale_height_r *
                         cello::kpc / units->length();
  double stellar_z   = enzo_config->method_background_acceleration_stellar_scale_height_z *
                         cello::kpc / units->length();
  double stellar_mass = enzo_config->method_background_acceleration_stellar_mass *
                        cello::mass_solar / units->mass();
  double bulge_mass   = enzo_config->method_background_acceleration_bulge_mass *
                        cello::mass_solar / units->mass();
  double bulgeradius = enzo_config->method_background_acceleration_bulge_radius *
                        cello::kpc / units->length();
  const double * amom = enzo_config->method_background_acceleration_angular_momentum;

  double G = this->G_four_pi_ *
             units->density() * units->time() * units->time();

  double rcore = enzo_config->method_background_acceleration_core_radius *
                 cello::kpc / units->length();

  //
  if (DM_mass > 0.0){
    double xtemp = DM_mass_radius / rcore;
    DM_density = DM_mass / ( (2.0 * cello::pi * pow(DM_mass_radius,3.0)) *
                              (0.5 * log(1.0 * xtemp * xtemp) + log(1.0 + xtemp) - atan(xtemp)));
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
         //double
         accel_sph = G * bulge_mass / pow(radius + bulgeradius,2) +    // bulge
                            G * DM_density * pow(rcore,3) *
                            (log(1.0+xtemp) - xtemp / (1.0+xtemp)) /
                            (radius * radius); // NFW DM profile

//                            cello::pi * G * DM_density * pow(rcore,3) / pow(radius,2)*
//                            (-2.0*atan(radius/rcore) + 2.0*log(1.0+radius/rcore)) +
//                               log(1.0+pow(radius/rcore,2));
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
         accel_R   = (rcyl    == 0.0 ? 0.0 : std::fabs(accel_R)   / (radius*cosmo_a));
         accel_z   = (zheight == 0.0 ? 0.0 : std::fabs(accel_z)*zheight/std::fabs(zheight) / cosmo_a);

         // now apply accelerations in cartesian (grid) coordinates
         int i = INDEX(ix,iy,iz,mx_,my_);
         if (ax) ax[i] -= (accel_sph * x + accel_R*xplane + accel_z*amom[0]);
         if (ay) ay[i] -= (accel_sph * y + accel_R*yplane + accel_z*amom[1]);
         if (az) az[i] -= (accel_sph * z + accel_R*zplane + accel_z*amom[2]);
        }
     }
  } // end loop over grid cells

  // Handle particles here - save for later

  return;
}
/*
double EnzoMethodBackgroundAcceleration::timestep (Block * block) const throw()
{
  return timestep_(block);
}
*/
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

  EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology * )
    block->simulation()->problem()->physics("cosmology");
  if (cosmology) {
    const int rank = block->rank();
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
