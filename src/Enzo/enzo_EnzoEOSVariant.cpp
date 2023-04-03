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

//----------------------------------------------------------------------

namespace {

std::vector<std::string> unexpected_keys_
(const std::unordered_map<std::string, double>& map,
 const std::vector<std::string>& expected) noexcept
{

  auto unexpected = [&expected](const std::string& s) -> bool
  {
    for (const std::string& e : expected) {if (s == e) {return false;}}
    return true;
  };

  std::vector<std::string> out;
  for (auto& kv : map) {if (unexpected(kv.first)) {out.push_back(kv.first);}}
  return out;
}

std::string strvec_to_string(const std::vector<std::string>& vec) noexcept
{
  std::string out = "[";
  for (std::size_t i = 0; i < vec.size(); i++) {
    if (i == 0) { out += ", "; }
    out += ("\"" + vec[i] + "\"");
  }
  return out + "]";
}

}

// to do: remove the following...
EnzoEOSVariant EnzoEOSVariant::factory
(const std::string& eos_type,
 const std::unordered_map<std::string, double>& params)
{
  if (eos_type == EnzoEOSIdeal::name()) {
    auto rslt = params.find("gamma");
    if (rslt == params.end()) {
      ERROR1("EnzoEOSVariant::factory",
             "eos type \"%s\" requires the \"gamma\" parameter",
             EnzoEOSIsothermal::name());
    } else if (params.size() != 1) {
      std::string tmp = strvec_to_string(unexpected_keys_(params,
                                                          {"gamma"}));
      ERROR2("EnzoEOSVariant::factory",
             "eos type \"%s\" doesn't expect the parameters: %s",
             EnzoEOSIsothermal::name(), tmp.c_str());
    }
    EnzoEOSIdeal tmp = {rslt->second};
    return EnzoEOSVariant(tmp);
  } else if (eos_type == EnzoEOSIsothermal::name()) {
    if (params.size() > 0) {
      ERROR1("EnzoEOSVariant::factory",
             "eos type \"%s\" doesn't expect any parameters",
             EnzoEOSIsothermal::name());
    }
    EnzoEOSIsothermal tmp{};
    return EnzoEOSVariant(tmp);
  } else {
    ERROR1("EnzoEOSVariant::factory",
           "No support is currently implemented for eos of the "
           "\"%s\" type",
           eos_type.c_str());
  }
}
