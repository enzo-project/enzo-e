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

EnzoMethodAccretionCompute::EnzoMethodAccretionCompute(double accretion_radius_cells)
  : Method(),
    accretion_radius_cells_(accretion_radius_cells)
{
  // This method requires three dimensions.
  ASSERT("EnzoMethodMergeStars::EnzoMethodMergeStars()",
	 "EnzoMethodMergeStars requires that we run a 3D problem (Domain: rank = 3)",
	 cello::rank());

  const EnzoConfig * enzo_config = enzo::config();

  const int * ghost_depth = enzo_config->field_ghost_depth;
  const int min_ghost_depth = std::min(ghost_depth[0],
				       std::min(ghost_depth[1],ghost_depth[2]));

  // The number of cell widths within the accretion radius cannot be larger than
  // any of the ghost depths
  ASSERT("EnzoMethodAccretionCompute::EnzoMethodAccretionCompute() ",
	 "The accretion radius must be no greater than the ghost depth"
	 "(4 by default)",
         accretion_radius_cells_ <= min_ghost_depth);

  // Check if merge_stars method is also being used and if it comes just before
  // accretion_compute in the method list
  const EnzoProblem * enzo_problem = enzo::problem();
  size_t i = 0;
  size_t merge_stars_index, accretion_compute_index;
  bool merge_stars_found = false;
  bool accretion_compute_found = false;
  
  while( enzo_problem->method(i) ){
    if ( enzo_problem->method(i)->name() == "merge_stars" )
      merge_stars_index = i;
    if ( enzo_problem->method(i)->name() == "accretion_compute" )
      accretion_compute_index = i;
    if (merge_stars_found && accretion_compute_found) break;
    ++i;
  }

  ASSERT("EnzoMethodAccretionCompute::EnzoMethodAccretionCompute() ",
	 "Accretion method requires running with MergeStars method also.",
         merge_stars_found);

  ASSERT("EnzoMethodAccretionCompute::EnzoMethodAccretionCompute() ",
	 "Accretion method requires running with MergeStars method also.",
	 accretion_compute_index == merge_stars_index + 1);

  // Check if merging radius is at least twice that of the accretion
  // radius
  ASSERT("EnzoMethodAccretionCompute::EnzoMethodAccretionCompute() ",
	 "Merging radius (Method:merge_stars:merging_radius_cells "
	 "must be at least twice the accretion radius "
	 "(Method:accretion_compute:accretion_radius).",
	 enzo_config->method_merge_stars_merging_radius_cells >=
	 2.0 * accretion_radius_cells_);

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
//   This does nothing at the moment - business is done in derived
//   classes

void EnzoMethodAccretionCompute::compute ( Block *block) throw()
{

  if (! block->is_leaf()) return;

  block->compute_done();

  return;
}

// Required
double EnzoMethodAccretionCompute::timestep ( Block *block) const throw()
{
  return std::numeric_limits<double>::max();
}
