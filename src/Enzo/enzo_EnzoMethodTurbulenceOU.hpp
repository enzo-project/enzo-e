// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulenceOU.hpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Wed Jul 23 00:31:13 UTC 2014
/// @date     Fri Aug 24 00:31:13 UTC 2018
/// @brief    [\ref Enzo] Implementation of Enzo IsoThermal TURBULENCE MHD method

#ifndef ENZO_ENZO_METHOD_TURBULENCE_OU_HPP
#define ENZO_ENZO_METHOD_TURBULENCE_OU_HPP

//----------------------------------------------------------------------

class EnzoMethodTurbulenceOU : public Method {

  /// @class    EnzoMethodTurbulenceOU
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ISOTHERMAL TURBULENCE MHD method

public: // interface

  /// Create a new EnzoMethodTurbulence object
    EnzoMethodTurbulenceOU
    (double gamma,
     const double domain_lower[3],
     const double domain_upper[3],
     bool apply_cooling,
     bool apply_forcing,
     bool apply_injection_rate,
     int cooling_term,
     double hc_alpha, 
     double hc_sigma,
     double injection_rate,
     double kfi,
     double kfa,
     double mach,
     int olap,
     bool read_sol,
     double sol_weight,
     double totemp,
     bool update_solution);

  /// Create an uninitialized EnzoMethodTurbulence object
  EnzoMethodTurbulenceOU()
    : Method()
  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodTurbulenceOU);

  /// Charm++ PUP::able migration constructor
  EnzoMethodTurbulenceOU (CkMigrateMessage *m)
    : Method (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  void compute_shift(EnzoBlock *, CkReductionMsg *msg);
  void compute_update(EnzoBlock *, CkReductionMsg *msg);

  virtual std::string name () throw () 
  { return "turbulence_ou"; }

private: // methods

private: // attributes

  double gamma_;
  bool apply_cooling_;
  bool apply_forcing_;
  bool apply_injection_rate_;
  int cooling_term_;
  double hc_alpha_; 
  double hc_sigma_;
  double injection_rate_;
  double kfi_;
  double kfa_;
  double mach_;
  int olap_;
  bool read_sol_;
  double sol_weight_;
  double totemp_;
  bool update_solution_;

  int is_NModes_;
  int is_Ndims_;

  // True only on first block's call on a node each cycle to avoid
  // excessive calls
  static int iupdate_phases_;

};

#endif /* ENZO_ENZO_METHOD_TURBULENCE_OU_HPP */
