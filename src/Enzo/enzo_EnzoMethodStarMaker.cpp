/// See LICENSE_CELLO file for license and copyright information

///
///
///
///

#include "cello.hpp"
#include "enzo.hpp"

//-------------------------------------------------------------------

EnzoMethodStarMaker::EnzoMethodStarMaker
()
  : Method()
{

  FieldDescr * field_descr = cello::field_descr();
  const EnzoConfig * enzo_config = enzo::config();

  // Initialize default Refresh object
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
                             enzo_sync_id_method_star_maker);
  refresh(ir)->add_all_fields();

  // Keep parametrers that will likely be relevant to all star
  // formation routines here
  use_density_threshold_ = enzo_config->method_star_maker_use_density_threshold;
  use_velocity_divergence_ = enzo_config->method_star_maker_use_velocity_divergence;
  use_dynamical_time_      = enzo_config->method_star_maker_use_dynamical_time;
  number_density_threshold_ = enzo_config->method_star_maker_number_density_threshold;
  efficiency_ = enzo_config->method_star_maker_efficiency;
  maximum_star_fraction_ = enzo_config->method_star_maker_maximum_mass_fraction;
  star_particle_mass_ = enzo_config->method_star_maker_minimum_star_mass;
}

//-------------------------------------------------------------------

void EnzoMethodStarMaker::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | use_density_threshold_;
  p | use_velocity_divergence_;
  p | use_dynamical_time_;
  p | number_density_threshold_;
  p | efficiency_;
  p | maximum_star_fraction_;
  p | star_particle_mass_;

}

//------------------------------------------------------------------
// REMOVE THIS COMPUTE ---- THIS IS JUST THE CLASS TEMPLATE
void EnzoMethodStarMaker::compute ( Block *block) throw()
{
  EnzoBlock * enzo_block = enzo::block(block);

  enzo_block->compute_done();

  return;
}

double EnzoMethodStarMaker::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodStarMaker::rescale_densities(EnzoBlock * enzo_block,
                                            const int index,
                                            const double density_ratio) throw() {

  // Loop through all passive scalars (colour fields)
  // which are mass fractions stored as densities, and rescale
  // to the new density after star formation
  // 
  //    density_ratio = old_density / new_density
  //

  Field field = enzo_block->data()->field();

  Grouping * field_groups = field.groups();
  int nc = field_groups->size("colour");

  for (int ic = 0; ic < nc; ic++){
    enzo_float * cfield = (enzo_float *)
      field.values(field_groups->item("colour",ic));

    cfield[index] *= density_ratio;

  }

  return;
}

/*
void EnzoMethodStarMaker::convert_densities_to_fraction(EnzoBlock * enzo_block,
                                                        int direction) throw() {

  // Actually, I don't really think we need this...
  //   this only needs to be done with cells that either have
  //   star formation or get gas removed for star formation. This is
  //   likely to be a small number of cells on a given grid, so
  //   really there is no need to do this conversion for every cell...


  Field field = enzo_block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0, &gx, &gy, &gz);

  int mx, my, mz;
  field.dimensions (0, &mx, &my, &mz);

  int nx, ny, nz;
  field.size ( &nx, &ny, &nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  enzo_float * d = (enzo_float *) field.values("density");

  if (direction == 0){ // convert density to fraction

    for (int iz = 0; iz < ngz; iz++){
      for (int iy = 0; iy < ngy; iy++){
        for (int ix = 0; ix < ngx; ix++){
          int i = INDEX(ix,iy,iz,ngx,ngy);

          double inv_dens = 1.0 / d[i];

          if (metal) metal[i] = metal[i] * inv_dens;
        }
      }
    }


  } else { // convert fraction to density

    if (metal) metal[i] *= d[i];
  }

  return;
}
*/

// ---------------------------------------------------------

int EnzoMethodStarMaker::check_number_density_threshold(
                                                       const double &d
                                                        ){

  ///  Apply the criteria that the local number density be greater
  ///  than the provided number density if use_density_threshold_ is
  ///  desired by the user.

  return !(use_density_threshold_) +
          (d >= number_density_threshold_);
}

int EnzoMethodStarMaker::check_velocity_divergence(
                double *vx, double *vy, double *vz,
                const int &index, const int &dix, const int &diy,
                const int &diz){

    ///  Apply the criteria that the divergence of the velocity
    ///  be negative, if so desired by user (use_velocity_divergence).

    if (!(this->use_velocity_divergence_)){
      return 1.0;
    }

    int result = 0;

    if(vx){
      result = (vx[index+dix] - vx[index-dix] < 0) ? 1 : 0;
    }

    if(vy && result){
      result = (vy[index+diy] - vy[index-diy] < 0) ? 1 : 0;
    }

    if(vz && result){
      result = (vz[index+diz] - vz[index-diz] < 0) ? 1 : 0;
    }

    return result;
}

int EnzoMethodStarMaker::check_minimum_mass(const double & m){
  /// Apply the condition that the mass of gas converted into
  /// stars in a single cell cannot exceed a certain fraction
  /// of that cell's mass
  return (maximum_star_fraction_ * m) > star_particle_mass_;
}
