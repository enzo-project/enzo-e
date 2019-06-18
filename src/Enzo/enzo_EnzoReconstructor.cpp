// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructor.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implements the EnzoReconstructor interface

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoReconstructor* EnzoReconstructor::construct_reconstructor
(std::vector<std::string> reconstructable_groups,
 std::vector<std::string> passive_groups, std::string name)
{
  // some repeated code from construct_riemann
  std::vector<std::string> groups = reconstructable_groups;
  groups.insert(groups.end(), passive_groups.begin(), passive_groups.end());

  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(name.size(), ' ');
  std::transform(name.begin(), name.end(), formatted.begin(),
		 ::tolower);
  EnzoReconstructor* out;
  if (formatted == std::string("nn")){
    out = new EnzoReconstructorNN(groups);
  } else if (formatted == std::string("plm")){
    out = new EnzoReconstructorPLM(groups);
  } else {
    ASSERT("EnzoReconstructor",
	   "The only allowed solvers are NN & PLM", false);
    out = NULL; // Deals with compiler warning
  }
  return out;
}
