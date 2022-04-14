/// See LICENSE_CELLO file for license and copyright information

/// @file   enzo_EnzoMethodAccretionCompute.cpp
/// @author Stefan Arridge (stefan.arridge@gmail.com)
/// @date   24 February 2022
/// @brief  Implementation of EnzoMethodAccretionCompute, a base class
///         for "accretion compute" methods. These methods compute
///         the accretion rate onto star particles, and change the properties
///         of the particles accordingly. Gas density is reduced by setting
///         negative values for the "density_accreted" field. The
///         "accretion_remove_gas" method then subtracts off density_accreted
///         from the gas density field


#include "cello.hpp"
#include "enzo.hpp"

EnzoMethodAccretionCompute::EnzoMethodAccretionCompute
(double accretion_radius_cells,
 double density_threshold)
  : Method(),
    accretion_radius_cells_(accretion_radius_cells),
    density_threshold_(density_threshold)
{
  // Check if density threshold is at least as large as the density floor
  // set by the VL+CT method
  ASSERT("EnzoMethodAccretionComputeDensThresh::EnzoMethodAccretionComputeDensThresh",
	 "Density threshold must be at least as large as the density "
	 "floor set by the VL+CT method",
         density_threshold_ >= enzo::config()->method_vlct_density_floor);

  // This method requires three dimensions.
  ASSERT("EnzoMethodMergeStars::EnzoMethodMergeStars()",
	 "EnzoMethodMergeStars requires that we run a 3D problem (Domain: rank = 3)",
	 cello::rank());

  ASSERT(
      "EnzoMethodAccretionCompute::EnzoMethodAccretionCompute()",
      "EnzoMethodAccretionCompute requires unigrid mode (Adapt : max_level = 0). "
      "In future, we may put in a refinement condition that blocks containing "
      "star particles are at the highest refinement level.",
      enzo::config()->mesh_max_level == 0);

  const int * ghost_depth = enzo::config()->field_ghost_depth;
  const int min_ghost_depth = std::min(ghost_depth[0],
				       std::min(ghost_depth[1],ghost_depth[2]));

  // The number of cell widths within the accretion radius cannot be larger than
  // any of the ghost depths
  ASSERT("EnzoMethodAccretionCompute::EnzoMethodAccretionCompute() ",
	 "The accretion radius must be no greater than the ghost depth"
	 "(4 by default)",
         accretion_radius_cells_ <= min_ghost_depth);

  // Refresh all fields
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();
}

void EnzoMethodAccretionCompute::pup (PUP::er &p)
{
  // NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | accretion_radius_cells_;
 
  return;
}


//------------------------------------------------------------------
//   This does nothing - business is done in derived
//   classes
void EnzoMethodAccretionCompute::compute(Block *block) throw() {

  block->compute_done();

  return;
}

// Required
double EnzoMethodAccretionCompute::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}

void EnzoMethodAccretionCompute::do_checks_() throw()
{
    // Check if merge_stars method precedes accretion_compute method
    ASSERT("EnzoMethodAccretionCompute",
	   "merge_stars must precede accretion_compute",
	   enzo::problem()->method_precedes("merge_stars",
					    "accretion_compute"));

    // Check if accretion_compute method precedes accretion_remove_gas
    // method
    ASSERT("EnzoMethodAccretionCompute",
	   "accretion_compute must precede accretion_remove_gas",
	   enzo::problem()->method_precedes("accretion_compute",
					    "accretion_remove_gas"));

    // Check if merging radius is at least twice that of the accretion
    // radius
    ASSERT("EnzoMethodAccretionCompute::EnzoMethodAccretionCompute() ",
	   "Merging radius (Method:merge_stars:merging_radius_cells "
	   "must be at least twice the accretion radius "
	   "(Method:accretion_compute:accretion_radius).",
	   enzo::config()->method_merge_stars_merging_radius_cells >=
	   2.0 * accretion_radius_cells_);

    // Check if VL+CT method is being used.
    ASSERT("EnzoMethodAccretionCompute::EnzoMethodAccretionCompute() ",
	   "accretion_compute requires vlct method. (note: at the "
	   "moment this means that cosmology can't be used in "
	   "combination with accretion.",
	   enzo::problem()->method_exists("vlct"));
}
