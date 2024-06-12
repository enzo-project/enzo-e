// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoPhysicsFluidProps.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2024-05-19
/// @brief    [\ref Enzo] Implementation of the EnzoPhysicsGravity class

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/gravity/gravity.hpp"
#include "parameters_ParameterGroup.hpp"

//----------------------------------------------------------------------

namespace{ // stuff inside an anonymous namespace is local to this file)

double parse_grav_const_codeU_(ParameterGroup p) {
  // for anybody looking at this implementation as an example of "How to access
  // parameters of a physics object," THIS IS A POOR EXAMPLE!!!
  // -> As a general rule, you should ONLY be accessing parameters within the
  //    associated parameter group (i.e. parameters directly accessible using
  //    the ``ParameterGroup`` argument). You should generally NOT use the
  //    global ``Parameters`` instance.
  // -> here we break these rules purely for the sake of maintaining backwards
  //    compatability!

  Parameters * all_parameters = cello::simulation()->parameters();

  std::string legacy_parname = "Method:gravity:grav_const";
  std::string actual_par_basename = "grav_const_codeU";

  const std::vector<std::string>& method_list = enzo::config()->method_list;
  bool has_grav_method = std::find(method_list.begin(), method_list.end(),
                            "gravity") != method_list.end();
  bool has_legacy_par = (has_grav_method &&
                         (all_parameters->param(legacy_parname) != nullptr));
  bool has_actual_par = p.param(actual_par_basename) != nullptr;

  if (has_legacy_par && has_actual_par) {
    ERROR3("parse_grav_const_codeU_",
           "\"%s\" isn't valid since \"%s:%s\" is specified.",
           legacy_parname.c_str(), p.get_group_path().c_str(),
           actual_par_basename.c_str());
  } else if (has_legacy_par) {
    WARNING3("parse_grav_const_codeU_",
             "\"%s\" is a legacy parameter that is replaced by \"%s:%s\"",
             legacy_parname.c_str(), p.get_group_path().c_str(),
             actual_par_basename.c_str());
    return all_parameters->value_float(legacy_parname, -1.0);
  } else {
    return p.value_float(actual_par_basename, -1.0);
  }

}

} // close anonymous namespace

//----------------------------------------------------------------------

EnzoPhysicsGravity::EnzoPhysicsGravity(ParameterGroup p)
  : Physics()
{

  // When grav_const_code_units is positive, it specifies the gravitational
  // constant in code units. Otherwise, the gravitational constant is
  // defined such that it's equal to the value of
  // enzo_constants::standard_grav_const when converted to cgs
  double grav_const_code_units = parse_grav_const_codeU_(p);

  if (grav_const_code_units <= 0.0) {
    grav_constant_codeU_ = -1.0;
  } else {
    ASSERT("EnzoPhysicsGravity::EnzoPhysicsGravity",
           "users aren't allowed to specify the gravitational constant in "
           "a cosmological simulation", enzo::cosmology() == nullptr);
    grav_constant_codeU_ = grav_const_code_units;
  }
}