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

namespace { // stuff inside anonymous namespace are local to this source file

struct BlockInfo {
  // In the future, we can make use of data::field_cells instead of this helper
  // struct (this helper struct only exists to avoid introducing changes to
  // test answers right now)


  // attributes:
  std::array<int, 3> dimensions;    // include ghost zones
  std::array<int, 3> ghost_depth;
  std::array<double, 3> cell_width;
  std::array<double, 3> lower_edge;

  BlockInfo(Block * block){
    Field field = block->data()->field();

    field.dimensions (0,&(dimensions[0]),&(dimensions[1]),&(dimensions[2]));
    field.ghost_depth(0,&(ghost_depth[0]),&(ghost_depth[1]),&(ghost_depth[2]));
    block->data()->lower(&(lower_edge[0]),&(lower_edge[1]),&(lower_edge[2]));
    double xp, yp, zp;
    block->data()->upper(&xp,&yp,&zp);
    field.cell_width(lower_edge[0],xp,&(cell_width[0]),
                     lower_edge[1],yp,&(cell_width[1]),
                     lower_edge[2],zp,&(cell_width[2]));
  }

  inline double x_val(int ix) const noexcept
  { return lower_edge[0] + (ix - ghost_depth[0] + 0.5) * cell_width[0]; }

  inline double y_val(int iy) const noexcept
  { return lower_edge[1] + (iy - ghost_depth[1] + 0.5) * cell_width[1]; }

  inline double z_val(int iz) const noexcept
  { return lower_edge[2] + (iz - ghost_depth[2] + 0.5) * cell_width[2]; }

};

//---------------------------------------------------------------------

struct GalaxyModelParameterPack {

  /// @class    GalaxyModelParameterPack
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] tracks the parameters needed for a background
  /// potential of a Galaxy model

  double DM_mass;
  double DM_mass_radius;
  double stellar_r;
  double stellar_z;
  double stellar_mass;
  double bulge_mass;
  double bulgeradius;
  std::array<double,3> amom;
  double rcore;

  static GalaxyModelParameterPack from_config(const EnzoConfig* enzo_config) {
    double DM_mass = enzo_config->method_background_acceleration_DM_mass;
    double DM_mass_radius = enzo_config->method_background_acceleration_DM_mass_radius;
    double stellar_r = enzo_config->method_background_acceleration_stellar_scale_height_r;
    double stellar_z = enzo_config->method_background_acceleration_stellar_scale_height_z;
    double stellar_mass = enzo_config->method_background_acceleration_stellar_mass;
    double bulge_mass = enzo_config->method_background_acceleration_bulge_mass;
    double bulgeradius = enzo_config->method_background_acceleration_bulge_radius;
    std::array<double,3> amom
      = {enzo_config->method_background_acceleration_angular_momentum[0],
         enzo_config->method_background_acceleration_angular_momentum[1],
         enzo_config->method_background_acceleration_angular_momentum[2]};
    double rcore = enzo_config->method_background_acceleration_core_radius;

    ASSERT1("GalaxyModel::GalaxyModel",
            "DM halo mass (=%e code_units) must be positive and specified in "
            "units of solar masses",
            DM_mass, (DM_mass > 0));

    return {DM_mass, DM_mass_radius, stellar_r, stellar_z,
            stellar_mass, bulge_mass, bulgeradius, amom, rcore};
  }

  void pup(PUP::er &p) {
    p | DM_mass;
    p | DM_mass_radius;
    p | stellar_r;
    p | stellar_z;
    p | stellar_mass;
    p | bulge_mass;
    p | bulgeradius;
    p | amom;
    p | rcore;
  }
};

class GalaxyModelFunctor {

  /// @class    GalaxyModelParameterPack
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Evaluates the background potential of a Galaxy model
  ///
  /// If we didn't support cosmological simulations, we would simply store the
  /// fields of the GalaxyModelParameterPack struct as members of this class
  /// (in terms of code units)

public:
  GalaxyModelFunctor(const GalaxyModelParameterPack& pack_dfltU,
                     const EnzoUnits* enzo_units)
  {
    pack_codeU_.DM_mass =
      pack_dfltU.DM_mass * enzo_constants::mass_solar / enzo_units->mass();
    pack_codeU_.DM_mass_radius =
      pack_dfltU.DM_mass_radius * enzo_constants::kpc_cm / enzo_units->length();
    pack_codeU_.stellar_r =
      pack_dfltU.stellar_r * enzo_constants::kpc_cm / enzo_units->length();
    pack_codeU_.stellar_z =
      pack_dfltU.stellar_z * enzo_constants::kpc_cm / enzo_units->length();
    pack_codeU_.stellar_mass =
      pack_dfltU.stellar_mass * enzo_constants::mass_solar / enzo_units->mass();
    pack_codeU_.bulge_mass =
      pack_dfltU.bulge_mass * enzo_constants::mass_solar / enzo_units->mass();
    pack_codeU_.bulgeradius =
      pack_dfltU.bulgeradius * enzo_constants::kpc_cm / enzo_units->length();
    pack_codeU_.amom = pack_dfltU.amom;
    pack_codeU_.rcore =
      pack_dfltU.rcore * enzo_constants::kpc_cm / enzo_units->length();

    double xtemp = pack_codeU_.DM_mass_radius / pack_codeU_.rcore;

    // compute the density constant for an NFW halo (rho_o)
    DM_density_ = (pack_codeU_.DM_mass / (4.0 * cello::pi * std::pow(pack_codeU_.rcore,3))) / (std::log(1.0+xtemp)-xtemp/(1.0+xtemp));
  }

  /// compute the x,y,z components of the acceleration
  ///
  /// @note
  /// it's unclear to me why we only consider the spherical acceleration
  /// component - this is equivalent to the way that the code was written
  /// before refactoring
  std::array<double,3> accel_fluid(double G_code, double cosmo_a,
                                   double x, double y, double z)
    const noexcept
  { return accel_helper_<true>(G_code, cosmo_a, x, y, z); }

  /// compute the x,y,z components of the acceleration
  std::array<double,3> accel_particle(double G_code, double cosmo_a,
                                      double x, double y, double z)
    const noexcept
  { return accel_helper_<false>(G_code, cosmo_a, x, y, z); }

private:

  template <bool only_sph_component>
  std::array<double,3> accel_helper_(double G_code, double cosmo_a,
                                     double x, double y, double z)
    const noexcept
  {
    double zheight = (pack_codeU_.amom[0]*x +
                      pack_codeU_.amom[1]*y +
                      pack_codeU_.amom[2]*z); // height above disk

    // projected positions in plane of the disk
    double xplane = x - zheight*pack_codeU_.amom[0];
    double yplane = y - zheight*pack_codeU_.amom[1];
    double zplane = z - zheight*pack_codeU_.amom[2];

    double radius = sqrt(xplane*xplane + yplane*yplane + zplane*zplane + zheight*zheight);
    double rcyl   = sqrt(xplane*xplane + yplane*yplane + zplane*zplane);

    // need to multiple all of the below by the gravitational constants
    double xtemp     = radius/pack_codeU_.rcore;

    double bulge_denom = pow(radius + pack_codeU_.bulgeradius, 2);
    double accel_sph =
      G_code * pack_codeU_.bulge_mass / bulge_denom +  // bulge
      4.0 * G_code * cello::pi * DM_density_ * pack_codeU_.rcore *
      (log(1.0+xtemp) - (xtemp / (1.0+xtemp))) / (xtemp*xtemp);

    double accel_R =
      G_code * pack_codeU_.stellar_mass * rcyl /
      sqrt( pow( pow(rcyl,2) +
                 pow(pack_codeU_.stellar_r +
                     sqrt( pow(zheight,2) + pow(pack_codeU_.stellar_z,2)),
                     2),
                 3)
            );

    double accel_z   =
      G_code * pack_codeU_.stellar_mass /
      sqrt(pow(zheight,2) + pow(pack_codeU_.stellar_z,2)) *
      zheight / sqrt(pow(pow(rcyl,2) +
                         pow(pack_codeU_.stellar_r +
                             sqrt(pow(zheight,2)
                                  + pow(pack_codeU_.stellar_z,2)),
                             2),
                         3))
      * (pack_codeU_.stellar_z * sqrt(pow(zheight,2) +
                                      pow(pack_codeU_.stellar_z,2)));

    accel_sph = (radius  == 0.0) ? 0.0 : std::fabs(accel_sph) / (radius*cosmo_a);
    accel_R   = (rcyl    == 0.0) ? 0.0 : std::fabs(accel_R)   / (rcyl*cosmo_a);
    accel_z   = (zheight == 0.0) ? 0.0 : std::fabs(accel_z)*zheight/std::fabs(zheight) / cosmo_a;

    if (only_sph_component) {
      accel_R = 0.0; accel_z = 0.0;
    }

    return { accel_sph * x + accel_R*xplane + accel_z*pack_codeU_.amom[0],  // x
             accel_sph * y + accel_R*yplane + accel_z*pack_codeU_.amom[1],  // y
             accel_sph * z + accel_R*zplane + accel_z*pack_codeU_.amom[2]}; // z
  }

private:
  GalaxyModelParameterPack pack_codeU_;

  /// density constant for an NFW halo (rho_o)
  double DM_density_;

};

//---------------------------------------------------------------------

void GalaxyModel(enzo_float * ax, enzo_float * ay, enzo_float * az,
                 double G_code, Particle * particle,
                 const BlockInfo block_info, const int rank,
                 const enzo_float cosmo_a,
                 const EnzoConfig * enzo_config,
                 const EnzoUnits * enzo_units, const double dt) throw();

//---------------------------------------------------------------------

void PointMass(enzo_float * ax, enzo_float * ay, enzo_float * az,
               double G_code, Particle * particle,
               const BlockInfo block_info, const int rank,
               const enzo_float cosmo_a,
               const EnzoConfig * enzo_config,
               const EnzoUnits * enzo_units, const double dt) throw();

} // close anonymous namespace

//---------------------------------------------------------------------

EnzoMethodBackgroundAcceleration::EnzoMethodBackgroundAcceleration
(bool zero_acceleration) // do need
 : Method(),
   zero_acceleration_(zero_acceleration)
{

  this->G_four_pi_ = 4.0 * cello::pi * enzo_constants::grav_constant;

  FieldDescr * field_descr = cello::field_descr();

  //const int id =  field_descr->field_id("density");
  const int iax = field_descr->field_id("acceleration_x");
  const int iay = field_descr->field_id("acceleration_y");
  const int iaz = field_descr->field_id("acceleration_z");


  // Do not need to refresh acceleration fields in this method
  // since we do not need to know any ghost zone information
  //
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field(iax);
  refresh->add_field(iay);
  refresh->add_field(iaz);

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

//---------------------------------------------------------------------

void EnzoMethodBackgroundAcceleration::compute_ (Block * block) throw()
{
  /* This applies background potential to grid and particles */

  //TRACE_METHOD("compute()",block);
  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();
  EnzoUnits * enzo_units = enzo::units();

  if (!(enzo_config->method_background_acceleration_apply_acceleration)) return;

  Field field = block->data()->field();

  BlockInfo block_info(block);

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
    if (rank >= 1) block_info.cell_width[0] *= cosmo_a;
    if (rank >= 2) block_info.cell_width[1] *= cosmo_a;
    if (rank >= 3) block_info.cell_width[2] *= cosmo_a;
  }


  enzo_float * ax = (enzo_float*) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float*) field.values ("acceleration_y");
  enzo_float * az = (enzo_float*) field.values ("acceleration_z");

  int m = (block_info.dimensions[0] * block_info.dimensions[1] *
           block_info.dimensions[2]);
  if (zero_acceleration_){
    if (ax){ for(int i = 0; i < m; i ++){ ax[i] = 0.0;}}
    if (ay){ for(int i = 0; i < m; i ++){ ay[i] = 0.0;}}
    if (az){ for(int i = 0; i < m; i ++){ az[i] = 0.0;}}
  }

  Particle particle = enzo_block->data()->particle();

  // unclear why the gravitational constant is different in each branch!

  if (enzo_config->method_background_acceleration_flavor == "GalaxyModel"){

    double G_code = enzo_constants::grav_constant * enzo_units->density() * enzo_units->time() * enzo_units->time();

    GalaxyModel(ax, ay, az, G_code, &particle, block_info, rank,
                cosmo_a, enzo_config, enzo_units, enzo_block->dt);

  } else if (enzo_config->method_background_acceleration_flavor == "PointMass"){

    double G_code = this->G_four_pi_ *
            enzo_units->density() * enzo_units->time() * enzo_units->time();
    PointMass(ax, ay, az, G_code, &particle, block_info, rank,
              cosmo_a, enzo_config, enzo_units, enzo_block->dt);
  } else {

    ERROR("EnzoMethodBackgroundAcceleration::compute_()",
          "Background acceleration flavor not recognized");

  }


  return;

} // compute

//---------------------------------------------------------------------

// define the main helper functions:
namespace { // stuff inside anonymous namespace are local to this source file

void PointMass(enzo_float * ax, enzo_float * ay, enzo_float * az,
               double G_code, Particle * particle,
               const BlockInfo block_info, const int rank,
               const enzo_float cosmo_a, const EnzoConfig * enzo_config,
               const EnzoUnits * enzo_units, const double dt) throw()
{
  // just need to define position of each cell

  double mass = enzo_config->method_background_acceleration_mass *
                enzo_constants::mass_solar / enzo_units->mass();
  double rcore = std::max(0.1*block_info.cell_width[0],
                     enzo_config->method_background_acceleration_core_radius/enzo_units->length());

  const double min_accel = mass / ((rcore*rcore*rcore)*cosmo_a);

  const int mx = block_info.dimensions[0];
  const int my = block_info.dimensions[1];
  const int mz = block_info.dimensions[2];

  const double* accel_center
    = enzo_config->method_background_acceleration_center;

  for (int iz=0; iz<mz; iz++){
    double z = (rank >= 3) ? block_info.z_val(iz) - accel_center[2] : 0.0;

    for (int iy=0; iy<my; iy++){
      double y = (rank >= 2) ? block_info.y_val(iy) - accel_center[1] : 0.0;

      for (int ix=0; ix<mx; ix++){
        double x = block_info.x_val(ix) - accel_center[0];

        double rsqr  = x*x + y*y + z*z;
        double r     = sqrt(rsqr);

        double accel = G_code * std::min(mass / ((rsqr)*r*cosmo_a), min_accel);

        int i = INDEX(ix,iy,iz,mx,my);

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

//---------------------------------------------------------------------

void GalaxyModel(enzo_float * ax, enzo_float * ay, enzo_float * az,
                 double G_code, Particle * particle,
                 const BlockInfo block_info, const int rank,
                 const enzo_float cosmo_a,
                 const EnzoConfig * enzo_config,
                 const EnzoUnits * enzo_units, const double dt) throw()
{
  double DM_mass     = enzo_config->method_background_acceleration_DM_mass *
                         enzo_constants::mass_solar / enzo_units->mass();
  double DM_mass_radius = enzo_config->method_background_acceleration_DM_mass_radius *
                          enzo_constants::kpc_cm / enzo_units->length();
  double stellar_r   = enzo_config->method_background_acceleration_stellar_scale_height_r *
                         enzo_constants::kpc_cm / enzo_units->length();
  double stellar_z   = enzo_config->method_background_acceleration_stellar_scale_height_z *
                         enzo_constants::kpc_cm / enzo_units->length();
  double stellar_mass = enzo_config->method_background_acceleration_stellar_mass *
                        enzo_constants::mass_solar / enzo_units->mass();
  double bulge_mass   = enzo_config->method_background_acceleration_bulge_mass *
                        enzo_constants::mass_solar / enzo_units->mass();
  double bulgeradius = enzo_config->method_background_acceleration_bulge_radius *
                        enzo_constants::kpc_cm / enzo_units->length();
  const double * amom = enzo_config->method_background_acceleration_angular_momentum;

  double rcore = enzo_config->method_background_acceleration_core_radius *
                 enzo_constants::kpc_cm / enzo_units->length();

  ASSERT1("Enzo::MethodBackgroundAcceleration", "DM halo mass (=%e code_units) must be positive and in units of solar masses", DM_mass, (DM_mass > 0));

  double xtemp = DM_mass_radius / rcore;

  // compute the density constant for an NFW halo (rho_o)
  double DM_density = (DM_mass / (4.0 * cello::pi * std::pow(rcore,3))) / (std::log(1.0+xtemp)-xtemp/(1.0+xtemp));

  double accel_sph, accel_R, accel_z;

  const int mx = block_info.dimensions[0];
  const int my = block_info.dimensions[1];
  const int mz = block_info.dimensions[2];

  const double* accel_center
    = enzo_config->method_background_acceleration_center;

  for (int iz=0; iz<mz; iz++){
     double z = (rank >= 3) ? block_info.z_val(iz) - accel_center[2] : 0.0;

     for (int iy=0; iy<my; iy++){
       double y = (rank >= 2) ? block_info.y_val(iy) - accel_center[1] : 0.0;

       for (int ix=0; ix<mx; ix++){
         double x = block_info.x_val(ix) - accel_center[0];

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
         accel_sph = G_code * bulge_mass / pow(radius + bulgeradius,2) +    // bulge
                     + 4.0 * G_code * cello::pi * DM_density * rcore *
                          (log(1.0+xtemp) - (xtemp / (1.0+xtemp))) / (xtemp*xtemp);

         //double
         accel_R   = G_code * stellar_mass * rcyl / sqrt( pow( pow(rcyl,2)
                             + pow(stellar_r + sqrt( pow(zheight,2)
                            + pow(stellar_z,2)),2),3));
         //double
         accel_z   = G_code * stellar_mass / sqrt(pow(zheight,2)
                               + pow(stellar_z,2))*zheight/sqrt(pow(pow(rcyl,2)
                               + pow(stellar_r + sqrt(pow(zheight,2)
                               + pow(stellar_z,2)),2),3))
                                   * (stellar_z * sqrt(pow(zheight,2) + pow(stellar_z,2)));

         accel_sph = (radius  == 0.0 ? 0.0 : std::fabs(accel_sph) / (radius*cosmo_a));
         accel_R   = (rcyl    == 0.0 ? 0.0 : std::fabs(accel_R)   / (rcyl*cosmo_a));
         accel_z   = (zheight == 0.0 ? 0.0 : std::fabs(accel_z)*zheight/std::fabs(zheight) / cosmo_a);

         accel_R = 0.0; accel_z = 0.0;
         // now apply accelerations in cartesian (grid) coordinates
         int i = INDEX(ix,iy,iz,mx,my);
         if (ax) ax[i] -= (accel_sph * x + accel_R*xplane + accel_z*amom[0]);
         if (ay) ay[i] -= (accel_sph * y + accel_R*yplane + accel_z*amom[1]);
         if (az) az[i] -= (accel_sph * z + accel_R*zplane + accel_z*amom[2]);
        }
     }
  } // end loop over grid cells

  // Update particle accelerations for gravitating particles

  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();

  int num_is_grav = particle_groups->size("is_gravitating");

  // Loop through particles to apply this to
  for (int ipt = 0; ipt < num_is_grav; ipt++){

    std::string particle_type = particle_groups->item("is_gravitating",ipt);
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

          const double x = px[ipdp] - accel_center[0];
          const double y = py[ipdp] - accel_center[1];
          const double z = pz[ipdp] - accel_center[2];

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
          accel_sph = G_code * bulge_mass / pow(radius + bulgeradius,2) +    // bulge
                     + 4.0 * G_code * cello::pi * DM_density * rcore *
                          (log(1.0+xtemp) - (xtemp / (1.0+xtemp))) / (xtemp*xtemp);

          accel_R   = G_code * stellar_mass * rcyl / sqrt( pow( pow(rcyl,2)
                              + pow(stellar_r + sqrt( pow(zheight,2)
                             + pow(stellar_z,2)),2),3));
          //double
          accel_z   = G_code * stellar_mass / sqrt(pow(zheight,2)
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

} // close anonymous namespace

//---------------------------------------------------------------------

double EnzoMethodBackgroundAcceleration::timestep (Block * block) throw()
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
