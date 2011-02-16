// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo PPM hydro method

#ifndef ENZO_ENZO_METHOD_PPM_HPP
#define ENZO_ENZO_METHOD_PPM_HPP

class EnzoMethodPpm : public Hyperbolic {

  /// @class    EnzoMethodPpm
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's PPM hydro method

public: // interface

  EnzoMethodPpm(Parameters * parameters,
		EnzoDescr  * enzo);

  /// Apply the method to advance a block one timestep 

  void compute_block(DataBlock * data_block,
		     double t, double dt) throw();

  /// Return the name of the method

  std::string name() const throw() 
  { return "ppm"; };

private:

  EnzoDescr * enzo_;

};

#endif /* ENZO_ENZO_METHOD_PPM_HPP */
