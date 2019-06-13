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

EnzoRiemann* EnzoRiemann::construct_riemann
(std::vector<std::string> integrable_groups,
 std::vector<std::string> passive_groups, std::string solver)
{
  // In the future, allocate array of flux functors here
  FluxFunctor** flux_funcs = NULL;
  int n_funcs = 0;

  // determine the type of solver to construct:
  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(solver.size(), ' ');
  std::transform(solver.begin(), solver.end(), formatted.begin(),
		 ::tolower);
  EnzoRiemann* out;

  // Eventually we may want do some check for non-MHD Riemann solvers

  if ((formatted == std::string("hlle")) ||
       (formatted == std::string("enzo_hlle"))){
    out = new EnzoRiemannHLLEEnzoMHD(integrable_groups, passive_groups,
				     flux_funcs, n_funcs);
  } else if (formatted == std::string("athena_hlle")){
    out = new EnzoRiemannHLLEAthenaMHD(integrable_groups, passive_groups,
				       flux_funcs, n_funcs);
  } else if (formatted == std::string("hlld")){
    // could possibly check that MHD fields are included
    out = new EnzoRiemannHLLD(integrable_groups, passive_groups, flux_funcs,
			      n_funcs);
  } else {
    ASSERT("EnzoRiemann",
	   "The only known solvers are HLLE (ENZO_HLLE), ATHENA_HLLE, & HLLD",
	   false);
    out = NULL;  // Deals with compiler warning
  }

  return out;
}
