// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSVariant.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-02-18
/// @brief    [\ref Enzo] Implementation of functions associated with the
///    EnzoEOSVariant class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

namespace { // stuff within anonymous namespace only visible in local file

struct NameVisitor{
  template<typename T> std::string operator()(T obj) { return T::name(); }
};

}

//----------------------------------------------------------------------

void pup(PUP::er &p, EnzoEOSVariant& variant) noexcept {

  // to make a hypothetical migration to std::variant as simple as possible, we
  // are explicitly avoiding the usage of macros

  const bool up = p.isUnpacking();

  // in C++ 14, we could write:
  // std::string eos_name = variant.visit( [](auto eos){return eos::name();} );
  std::string eos_name = variant.visit( NameVisitor() );

  p | eos_name;

  // We make the fundamental assumption that EOS objects of different types
  // return unique names

  if (EnzoEOSIdeal::name() == eos_name) {
    if (up) { variant = EnzoEOSVariant(EnzoEOSIdeal()); }
    pup(p, variant.get<EnzoEOSIdeal>());

  } else if (EnzoEOSIsothermal::name() == eos_name) {
    if (up) { variant = EnzoEOSVariant(EnzoEOSIsothermal()); }
    pup(p, variant.get<EnzoEOSIsothermal>());

  } else {
    ERROR1("pup(PUP::er&, EnzoEOSVariant&)",
           "There is currently no support for PUPing the EOS alternative with "
           "type \"%s\". Did you add a new EOS type and forget to update the "
           "PUP routine?",
           eos_name.c_str());
  }
}
