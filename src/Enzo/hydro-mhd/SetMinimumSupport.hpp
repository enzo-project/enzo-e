// See LICENSE_CELLO file for license and copyright information

/// @file     SetMinimumSupport.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-07-21
/// @brief    Declare function to set the energy to provide minimal pressure
///           support
///
/// When this function was first ported to Enzo-E, it was an instance method of
/// the EnzoBlock class. In the future, it probably makes sense to integrate
/// this function into the fluid-props functionality.
///
/// In the short-term, we leave it as a free-standing function

#ifndef ENZO_HYDROMHD_SETMINIMUMSUPPORT_HPP
#define ENZO_HYDROMHD_SETMINIMUMSUPPORT_HPP

namespace enzo{

  /// Set the energy to provide minimal pressure support
  int SetMinimumSupport(EnzoBlock& block,
                        enzo_float &MinimumSupportEnergyCoefficient,
                        enzo_float minimum_pressure_support_parameter,
                        bool comoving_coordinates);
}

#endif /* ENZO_HYDROMHD_SETMINIMUMSUPPORT_HPP */
