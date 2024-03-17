// See LICENSE_CELLO file for license and copyright information

/// @file     parameters_ParameterGroup.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Jan 17 2023
/// @brief    [\ref Parameters] Implementation for the ParameterGroup class

#include "cello.hpp"

#include "parameters.hpp"

//----------------------------------------------------------------------

Param * ParameterGroup::param (std::string parameter)
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  Param* out = wrapped_p_.param(full_name(parameter));
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

parameter_type ParameterGroup::type(std::string param) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  parameter_type out = wrapped_p_.type(full_name(param));
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

int ParameterGroup::value_integer (std::string s, int deflt) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  int out = wrapped_p_.value_integer(full_name(s), deflt);
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

double ParameterGroup::value_float (std::string s, double deflt) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  double out = wrapped_p_.value_float(full_name(s), deflt);
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

bool ParameterGroup::value_logical (std::string s, bool deflt) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  bool out = wrapped_p_.value_logical(full_name(s), deflt);
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

std::string ParameterGroup::value_string (std::string s,
                                          std::string deflt) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  std::string out = wrapped_p_.value_string(full_name(s), deflt);
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

int ParameterGroup::list_length (std::string parameter)
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  int out = wrapped_p_.list_length(full_name(parameter));
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

int ParameterGroup::list_value_integer (int i, std::string s,
                                        int deflt) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  int out = wrapped_p_.list_value_integer(i, full_name(s), deflt);
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

double ParameterGroup::list_value_float (int i, std::string s,
                                         double deflt) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  double out = wrapped_p_.list_value_float(i, full_name(s), deflt);
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

bool ParameterGroup::list_value_logical (int i, std::string s,
                                         bool deflt) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  bool out = wrapped_p_.list_value_logical(i, full_name(s), deflt);
  restore_wrapped_p_groups_(init_groups);
  return out;
}

//----------------------------------------------------------------------

std::string ParameterGroup::list_value_string (int i, std::string s,
                                               std::string deflt) noexcept
{
  std::vector<std::string> init_groups(pop_wrapped_p_groups_());
  std::string out = wrapped_p_.list_value_string(i, full_name(s), deflt);
  restore_wrapped_p_groups_(init_groups);
  return out;
}
