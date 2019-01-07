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

    int result = 0.0;

    if(vx){
      result += vx[index+dix] - vx[index-dix];
    }

    if(vy){
      result += vy[index+diy] - vy[index-diy];
    }

    if(vz){
      result += vz[index+diz] - vz[index-diz];
    }

    return result;
}

int EnzoMethodStarMaker::check_minimum_mass(const double & m){
  /// Apply the condition that the mass of gas converted into
  /// stars in a single cell cannot exceed a certain fraction
  /// of that cell's mass
  return (maximum_star_fraction_ * m) > star_particle_mass_;
}
