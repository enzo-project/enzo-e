// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructor.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implements the EnzoReconstructor interface

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoReconstructor* EnzoReconstructor::construct_reconstructor
(const std::vector<std::string> &active_reconstructed_quantities,
 std::string name, enzo_float theta_limiter)
{
  
  ASSERT("EnzoReconstructor::construct_reconstructor",
	 "theta_limiter must satisfy 1<=theta_limiter<=2",
	 (1.<=theta_limiter) && (theta_limiter<=2));

  // Construct a vector of keys for all components of the quantities listed in
  // active_reconstructed_quantities
  std::vector<std::string> key_l;

  for (const std::string &quantity_name : active_reconstructed_quantities){
    bool vector_quantity;
    if (!EnzoCenteredFieldRegistry::quantity_properties(quantity_name,
                                                        &vector_quantity)){
      ERROR1("EnzoReconstructorNN::reconstruct_interface",
             "\"%s\" is not a known quantity.", quantity_name.c_str());
    }
    if (vector_quantity){
      key_l.push_back(quantity_name + "_x");
      key_l.push_back(quantity_name + "_y");
      key_l.push_back(quantity_name + "_z");
    } else {
      key_l.push_back(quantity_name);
    }
  }

  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(name.size(), ' ');
  std::transform(name.begin(), name.end(), formatted.begin(),
		 ::tolower);
  EnzoReconstructor* out;
  if (formatted == std::string("nn")){
    out = new EnzoReconstructorNN(key_l);
  } else if ((formatted == std::string("plm")) ||
	     (formatted == std::string("plm_enzo"))){
    out = new EnzoReconstructorPLMEnzoRKLim(key_l, theta_limiter);
  } else if (formatted == std::string("plm_athena")) {
    out = new EnzoReconstructorPLMAthenaLim(key_l, theta_limiter);
  } else {
    ASSERT("EnzoReconstructor",
	   "The only allowed solvers are NN, PLM, PLM_ENZO, & PLM_ATHENA",
	   false);
    out = NULL; // Deals with compiler warning
  }
  return out;
}
