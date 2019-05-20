// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemann.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of EnzoRiemann

#include <string>
#include <algorithm>
#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoRiemann* EnzoRiemann::construct_riemann(std::string solver,
					    const EnzoFieldConditions cond)
{
  // Generate vector of group names of passive advected scalars
  EnzoCenteredFieldRegistry registry;
  std::vector<std::string> passive_groups;
  passive_groups = registry.passive_scalar_group_names();

  // In the future, allocate array of flux functors here
  FluxFunctor** flux_funcs = NULL;
  int n_funcs = 0;

  // determine the type of solver to construct:
  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(solver.size(), ' ');
  std::transform(solver.begin(), solver.end(), formatted.begin(),
		 ::tolower);
  EnzoRiemann* out;

  if (formatted == std::string("hlle")){
    out = new EnzoRiemannHLLE(cond, passive_groups, flux_funcs, n_funcs);
  } else if (formatted == std::string("hlld")){
    out = new EnzoRiemannHLLD(cond, passive_groups, flux_funcs, n_funcs);
  } else {
    ASSERT("EnzoRiemann", "The only allowed solvers are HLLE & HLLD", false);
    out = NULL;  // Deals with compiler warning
  }

  return out;
}
