// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemann.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of EnzoRiemann

#include <string>
#include <algorithm>
#include "cello.hpp"
#include "enzo.hpp"

// public header:
#include "EnzoRiemann.hpp"

// private headers:
#include "EnzoRiemannLUT.hpp"
#include "EnzoRiemannUtils.hpp"
#include "EnzoRiemannImpl.hpp"
#include "EnzoRiemannHLL.hpp"
#include "EnzoRiemannHLLC.hpp"
#include "EnzoRiemannHLLD.hpp"

//----------------------------------------------------------------------

EnzoRiemann* EnzoRiemann::construct_riemann
(const EnzoRiemann::FactoryArgs& factory_args) noexcept
{
  // determine the type of solver to construct:
  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(factory_args.solver.size(), ' ');
  std::transform(factory_args.solver.begin(), factory_args.solver.end(),
                 formatted.begin(), ::tolower);
  EnzoRiemann* out = nullptr; // set to NULL to suppress compiler warnings

  if (formatted == "hll"){
    ASSERT("EnzoRiemann::construct_riemann",
           ("An \"HLL\" Riemann solver without magnetic fields isn't "
            "currently supported."), factory_args.mhd);
    out = new EnzoRiemannHLLMHD(factory_args, factory_args.internal_energy);

  } else if (formatted == "hlle"){
    if (factory_args.mhd){
      out = new EnzoRiemannHLLEMHD(factory_args, factory_args.internal_energy);
    } else {
      ERROR("EnzoRiemann::construct_riemann",
            "The \"HLLE\" Riemann solver without magnetic fields is untested");
      out = new EnzoRiemannHLLE(factory_args, factory_args.internal_energy);
    }

  } else if (formatted == std::string("hllc")){
    ASSERT("EnzoRiemann::construct_riemann",
           "The \"HLLC\" Riemann Solver can't support mhd", !factory_args.mhd);
    out = new EnzoRiemannHLLC(factory_args, factory_args.internal_energy);

  } else if (formatted == std::string("hlld")){
    ASSERT("EnzoRiemann::construct_riemann",
           "The \"HLLD\" Riemann Solver requires magnetic fields",
           factory_args.mhd);
    out = new EnzoRiemannHLLD(factory_args, factory_args.internal_energy);

  } else {
    ERROR("EnzoRiemann::construct_riemann",
	  "The only known solvers are HLL, HLLE, HLLC, & HLLD");
  }

  return out;
}

//----------------------------------------------------------------------

void EnzoRiemann::check_key_order_(const EnzoEFltArrayMap &map, bool prim,
                                   const str_vec_t &passive_list) const noexcept
{
  str_vec_t standard_keys =
    prim ? primitive_quantity_keys() : integration_quantity_keys();
  // confirm that the first `standard_keys.size()` keys of map have the same
  // order as `standard_keys`
  map.validate_key_order(standard_keys, true, true);

  // confirm that all keys in passive_list also appear in map
  // (for now, we'll be permissive about the order of these keys)
  for (const auto& key : passive_list){
    if (!map.contains(key)){
      const std::string &name = map.name();
      ERROR2("EnzoRiemannImpl<ImplFunctor>::check_key_order_",
             "The \"%s\" map is missing a \"%s\" passive scalar key",
             name.c_str(), key.c_str());
    }
  }
}
