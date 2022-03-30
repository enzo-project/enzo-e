// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructor.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implements the EnzoReconstructor interface

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoReconstructor* EnzoReconstructor::construct_reconstructor
(const std::vector<std::string> &active_primitive_keys,
 std::string name, enzo_float theta_limiter)
{
  
  ASSERT("EnzoReconstructor::construct_reconstructor",
	 "theta_limiter must satisfy 1<=theta_limiter<=2",
	 (1.<=theta_limiter) && (theta_limiter<=2));

  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(name.size(), ' ');
  std::transform(name.begin(), name.end(), formatted.begin(),
		 ::tolower);
  EnzoReconstructor* out;
  if (formatted == std::string("nn")){
    out = new EnzoReconstructorNN(active_primitive_keys);
  } else if ((formatted == std::string("plm")) ||
	     (formatted == std::string("plm_enzo"))){
    out = new EnzoReconstructorPLMEnzoRKLim(active_primitive_keys,
					    theta_limiter);
  } else if (formatted == std::string("plm_athena")) {
    out = new EnzoReconstructorPLMAthenaLim(active_primitive_keys,
					    theta_limiter);
  } else {
    ASSERT("EnzoReconstructor",
	   "The only allowed solvers are NN, PLM, PLM_ENZO, & PLM_ATHENA",
	   false);
    out = NULL; // Deals with compiler warning
  }
  return out;
}
