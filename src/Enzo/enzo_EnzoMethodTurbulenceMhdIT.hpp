// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulenceMhdIT.hpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Wed Jul 23 00:31:13 UTC 2014
/// @date     Fri Aug 24 00:31:13 UTC 2018
/// @brief    [\ref Enzo] Implementation of Enzo IsoThermal TURBULENCE MHD method

#ifndef ENZO_ENZO_METHOD_TURBULENCE_MHD_IT_HPP
#define ENZO_ENZO_METHOD_TURBULENCE_MHD_IT_HPP

//----------------------------------------------------------------------

class EnzoMethodTurbulenceMhdIT : public Method {

  /// @class    EnzoMethodTurbulenceMhdIT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ISOTHERMAL TURBULENCE MHD method

public: // interface

  /// Create a new EnzoMethodTurbulence object
  EnzoMethodTurbulenceMhdIT(double edot,
			    double density_initial,
			    double bfieldx_initial,
			    double mach_number,
			    bool comoving_coordinates);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodTurbulenceMhdIT);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodTurbulenceMhdIT (CkMigrateMessage *m)
    : Method (m),
      density_initial_(0.0),
      bfieldx_initial_(0.0),
      edot_(0.0),
      mach_number_(0.0),
      comoving_coordinates_(false)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "turbulence_mhd_it"; }

  /// Resume computation after a reduction
  virtual void compute_resume ( Block * block,
				CkReductionMsg * msg) throw(); 

private: // methods

  void compute_resume_ (Block * block, CkReductionMsg * msg) throw();
  void monitor_output_(Block * block, double * g, double norm, double bnotx);

private: // attributes

  // Initial density
  double density_initial_;

  // Initial B-field
  double bfieldx_initial_;

  // Corresponds to Enzo "RandomForcingEdot" parameter
  double edot_;

  // Mach number
  double mach_number_;

  // Comoving Coordinates
  bool comoving_coordinates_;
};

#endif /* ENZO_ENZO_METHOD_TURBULENCE_MHD_IT_HPP */
