// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulence.hpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Wed Jul 23 00:31:13 UTC 2014
/// @brief    [\ref Enzo] Implementation of Enzo TURBULENCE hydro method

#ifndef ENZO_ENZO_METHOD_TURBULENCE_HPP
#define ENZO_ENZO_METHOD_TURBULENCE_HPP

//----------------------------------------------------------------------

enum {
  index_turbulence_vad,
  index_turbulence_aad,
  index_turbulence_vvdot,
  index_turbulence_vvot,
  index_turbulence_vvd,
  index_turbulence_vv,
  index_turbulence_dd,
  index_turbulence_d,
  index_turbulence_dax,
  index_turbulence_day,
  index_turbulence_daz,
  index_turbulence_dvx,
  index_turbulence_dvy,
  index_turbulence_dvz,
  index_turbulence_dlnd,
  index_turbulence_zones,
  index_turbulence_mind,
  index_turbulence_maxd,
  max_turbulence_array };

//----------------------------------------------------------------------

class EnzoMethodTurbulence : public Method {

  /// @class    EnzoMethodTurbulence
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's TURBULENCE hydro method

public: // interface

  /// Create a new EnzoMethodTurbulence object
  EnzoMethodTurbulence(double edot,
		       double density_initial,
		       double temperature_initial,
		       double mach_number,
		       bool comoving_coordinates);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodTurbulence);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodTurbulence (CkMigrateMessage *m)
    : Method (m),
      density_initial_(0.0),
      temperature_initial_(0.0),
      edot_(0.0),
      mach_number_(0.0),
      comoving_coordinates_(false)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "turbulence"; }

  /// Resume computation after a reduction
  virtual void compute_resume ( Block * block,
				CkReductionMsg * msg) throw(); 

private: // methods

  void compute_resume_ (Block * block, CkReductionMsg * msg) throw();

private: // attributes

  // Initial density
  double density_initial_;

  // Initial temperature
  double temperature_initial_;

  // Corresponds to Enzo "RandomForcingEdot" parameter
  double edot_;

  // Mach number
  double mach_number_;

  // Comoving Coordinates
  bool comoving_coordinates_;
};

#endif /* ENZO_ENZO_METHOD_TURBULENCE_HPP */
