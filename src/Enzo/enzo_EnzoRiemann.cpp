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

EnzoRiemann* EnzoRiemann::construct_riemann(const std::string& solver, const bool mhd,
                                            const bool internal_energy)
{
  // determine the type of solver to construct:
  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(solver.size(), ' ');
  std::transform(solver.begin(), solver.end(), formatted.begin(),
		 ::tolower);
  EnzoRiemann* out = nullptr; // set to NULL to suppress compiler warnings

  if (formatted == "hll"){
    ASSERT("EnzoRiemann::construct_riemann",
           ("An \"HLL\" Riemann solver without magnetic fields isn't "
            "currently supported."), mhd);
    out = new EnzoRiemannHLLMHD(internal_energy);

  } else if (formatted == "hlle"){
    if (mhd){
      out = new EnzoRiemannHLLEMHD(internal_energy);
    } else {
      ERROR("EnzoRiemann::construct_riemann",
            "The \"HLLE\" Riemann solver without magnetic fields is untested");
      out = new EnzoRiemannHLLE(internal_energy);
    }

  } else if (formatted == std::string("hllc")){
    ASSERT("EnzoRiemann::construct_riemann",
           "The \"HLLC\" Riemann Solver can't support mhd", !mhd);
    out = new EnzoRiemannHLLC(internal_energy);

  } else if (formatted == std::string("hlld")){
    ASSERT("EnzoRiemann::construct_riemann",
           "The \"HLLD\" Riemann Solver requires magnetic fields", mhd);
    out = new EnzoRiemannHLLD(internal_energy);

  } else {
    ERROR("EnzoRiemann::construct_riemann",
	  "The only known solvers are HLL, HLLE, HLLC, & HLLD");
  }

  return out;
}
