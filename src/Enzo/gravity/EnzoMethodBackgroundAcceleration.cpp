// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodBackgroundAcceleration.cpp
/// @author   Andrew Emerick (aemerick11@gmail.com)
/// @date     2018-05
/// @brief    Implements the EnzoMethodBackgroundAcceleration class

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/gravity/gravity.hpp"

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

class GalaxyModel {

  /// @class    GalaxyModel
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Evaluates the background potential of a Galaxy model
  ///
  /// If we didn't support cosmological simulations, we would simply store the
  /// fields of the EnzoPotentialConfigGalaxy struct as members of this class
  /// (in terms of code units)

public:
  GalaxyModel(const EnzoPotentialConfigGalaxy& pack_dfltU, const Units* units)
  {
    pack_codeU_ = EnzoPotentialConfigGalaxy::to_codeU(pack_dfltU, units);

    double xtemp = pack_codeU_.DM_mass_radius / pack_codeU_.rcore;

    // compute the density constant for an NFW halo (rho_o)
    DM_density_ = (pack_codeU_.DM_mass / (4.0 * cello::pi * std::pow(pack_codeU_.rcore,3))) / (std::log(1.0+xtemp)-xtemp/(1.0+xtemp));
  }

  /// compute the x,y,z components of the acceleration
  std::array<double,3> accel(double G_code, double cosmo_a,
                             double x, double y, double z)
    const noexcept
  {
    double zheight = (pack_codeU_.amom_uvec[0]*x +
                      pack_codeU_.amom_uvec[1]*y +
                      pack_codeU_.amom_uvec[2]*z); // height above disk

    // projected positions in plane of the disk
    double xplane = x - zheight*pack_codeU_.amom_uvec[0];
    double yplane = y - zheight*pack_codeU_.amom_uvec[1];
    double zplane = z - zheight*pack_codeU_.amom_uvec[2];

    double radius = sqrt(xplane*xplane + yplane*yplane + zplane*zplane + zheight*zheight);
    double rcyl   = sqrt(xplane*xplane + yplane*yplane + zplane*zplane);

    // need to multiple all of the below by the gravitational constants
    double xtemp     = radius/pack_codeU_.rcore;

    double bulge_denom = pow(radius + pack_codeU_.bulge_radius, 2);
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

    return {accel_sph * x + accel_R*xplane + accel_z*pack_codeU_.amom_uvec[0],
            accel_sph * y + accel_R*yplane + accel_z*pack_codeU_.amom_uvec[1],
            accel_sph * z + accel_R*zplane + accel_z*pack_codeU_.amom_uvec[2]};
  }

private:
  EnzoPotentialConfigGalaxy pack_codeU_;

  /// density constant for an NFW halo (rho_o)
  double DM_density_;

};

//---------------------------------------------------------------------

class PointMassModel {

public:

  PointMassModel(const EnzoPotentialConfigPointMass& pack_dfltU,
                 const Units* units, double cosmo_a,
                 std::array<double,3> cell_width)
  {
    pack_codeU_ = EnzoPotentialConfigPointMass::to_codeU(pack_dfltU, units);

    // TODO: should we be using the min or max cell_width?
    double rcore_tmp = std::max(0.1*cell_width[0], pack_codeU_.rcore);
    min_accel_ = pack_codeU_.mass / ((rcore_tmp*rcore_tmp*rcore_tmp)*cosmo_a);
  }

  std::array<double,3> accel(double G_code, double cosmo_a,
                             double x, double y, double z)
    const noexcept
  {
    double rsqr  = x*x + y*y + z*z;
    double r     = sqrt(rsqr);

    double accel =
      G_code * std::min(pack_codeU_.mass / ((rsqr)*r*cosmo_a), min_accel_);

    return {accel * x, accel * y, accel * z};
  }

private:
  EnzoPotentialConfigPointMass pack_codeU_;
  double min_accel_;
};

//---------------------------------------------------------------------


template<typename T>
void compute_accel_(const T functor,
                    enzo_float * ax, enzo_float * ay, enzo_float * az,
                    double G_code, Particle * particle,
                    const BlockInfo block_info, const int rank,
                    const enzo_float cosmo_a,
                    const std::array<double, 3> accel_center,
                    const double dt) noexcept
{
  const int mx = block_info.dimensions[0];
  const int my = block_info.dimensions[1];
  const int mz = block_info.dimensions[2];

  for (int iz=0; iz<mz; iz++){
     double z = (rank >= 3) ? block_info.z_val(iz) - accel_center[2] : 0.0;

     for (int iy=0; iy<my; iy++){
       double y = (rank >= 2) ? block_info.y_val(iy) - accel_center[1] : 0.0;

       for (int ix=0; ix<mx; ix++){
         double x = block_info.x_val(ix) - accel_center[0];

         std::array<double,3> accel = functor.accel(G_code, cosmo_a, x, y, z);
         // now apply accelerations in cartesian (grid) coordinates
         int i = INDEX(ix,iy,iz,mx,my);
         if (ax) ax[i] -= accel[0];
         if (ay) ay[i] -= accel[1];
         if (az) az[i] -= accel[2];
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

          std::array<double,3> accel = functor.accel(G_code, cosmo_a, x, y, z);

          // now apply accelerations in cartesian (grid) coordinates
          if (pax) pax[ipda] -= accel[0];
          if (pay) pay[ipda] -= accel[1];
          if (paz) paz[ipda] -= accel[2];

        } // end loop over particles

      } // end loop over batch

    } // end if particles exist

  } // end loop over particles

}

} // close anonymous namespace

//---------------------------------------------------------------------

EnzoMethodBackgroundAcceleration::EnzoMethodBackgroundAcceleration
(ParameterGroup p)
 : Method(),
   zero_acceleration_(false),
   potential_center_xyz_{}, // fills array with zeros
   flavor_(p.value_string("flavor","unknown")),
   galaxy_pack_dfltU_(nullptr),
   point_mass_pack_dfltU_(nullptr)
{

  // If self-gravity is calculated, we do not need to zero out the acceleration
  // field from the previous time stepbefore adding the background acceleration
  bool preceded_by_gravity = enzo::problem()->method("gravity") != nullptr;
  zero_acceleration_ = !preceded_by_gravity;

  for (int i = 0; i < 3; i++) {
    potential_center_xyz_[i] = p.list_value_float(i,"center",0.5);
  }

  if (flavor_ == "GalaxyModel") {
    galaxy_pack_dfltU_ = std::unique_ptr<EnzoPotentialConfigGalaxy>
      (new EnzoPotentialConfigGalaxy);
    *galaxy_pack_dfltU_ = EnzoPotentialConfigGalaxy::from_parameters(p);

  } else if (flavor_ == "PointMass") {
    point_mass_pack_dfltU_ = std::unique_ptr<EnzoPotentialConfigPointMass>
      (new EnzoPotentialConfigPointMass);
    *point_mass_pack_dfltU_ = EnzoPotentialConfigPointMass::from_parameters(p);

  } else {
    ERROR1("EnzoMethodBackgroundAcceleration::EnzoMethodBackgroundAcceleration",
           "Background acceleration flavor not recognized: %s",
           flavor_.c_str());
  }

  FieldDescr * field_descr = cello::field_descr();

  const int iax = field_descr->field_id("acceleration_x");
  const int iay = field_descr->field_id("acceleration_y");
  const int iaz = field_descr->field_id("acceleration_z");


  // Do not need to refresh acceleration fields in this method
  // since we do not need to know any ghost zone information
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field(iax);
  refresh->add_field(iay);
  refresh->add_field(iaz);

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
  Units * units = (Units*)enzo::units();

  Field field = block->data()->field();

  BlockInfo block_info(block);

  // We will probably never be in the situation of constant acceleration
  // and cosmology, but just in case.....
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  enzo_float cosmo_a = 1.0;

  const int rank = cello::rank();
  const double dt   = block->state()->dt();
  const double time = block->state()->time();

  if (cosmology) {
    enzo_float cosmo_dadt = 0.0;
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

  if (galaxy_pack_dfltU_ != nullptr) {

    double G_code = enzo::grav_constant_cgs() * units->density() * units->time() * units->time();
    const GalaxyModel functor(*galaxy_pack_dfltU_, units);

    compute_accel_(functor, ax, ay, az, G_code, &particle, block_info, rank,
                   cosmo_a, potential_center_xyz_, dt);

  } else if (point_mass_pack_dfltU_ != nullptr) {

    double G_code = (4.0 * cello::pi * enzo::grav_constant_cgs()) *
            units->density() * units->time() * units->time();
    const PointMassModel functor(*point_mass_pack_dfltU_, units,
                                 cosmo_a, block_info.cell_width);

    compute_accel_(functor, ax, ay, az, G_code, &particle, block_info, rank,
                   cosmo_a, potential_center_xyz_, dt);

  } else {

    ERROR("EnzoMethodBackgroundAcceleration::compute_()",
          "Background acceleration flavor not recognized");

  }


  return;

} // compute

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
    const double dt   = block->state()->dt();
    const double time = block->state()->time();
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
