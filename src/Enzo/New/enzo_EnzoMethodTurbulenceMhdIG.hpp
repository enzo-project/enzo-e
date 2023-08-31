// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMHDTurbulenceIG.hpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Wed Jul 23 00:31:13 UTC 2014
/// @date     Thu Sep 20 00:31:13 UTC 2018
/// @brief    [\ref Enzo] Implementation of Enzo Ideal Gas TURBULENCE MHD method with Ornstein-Uhlenbeck pumping

#ifndef ENZO_ENZO_METHOD_MHDTURBULENCEIG_HPP
#define ENZO_ENZO_METHOD_MHDTURBULENCEIG_HPP

//----------------------------------------------------------------------

class EnzoMethodTurbulenceMhdIG : public Method {

  /// @class    EnzoMethodTurbulenceMhdIG
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's IDEAL GAS MHD TURBULENCE forcing method

public: // interface

  /// Create a new EnzoMethodTurbulenceIG object
  EnzoMethodTurbulenceMhdIG
  (double gamma,
   double density_initial,
   double pressure_initial,
   double bfieldx_initial,
   double mach_number,
   double solenoidal_fraction,
   double kfmin,
   double kfmax,
   bool comoving_coordinates);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodTurbulenceMhdIG);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodTurbulenceMhdIG (CkMigrateMessage *m)
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
  { return "turbulence_mhd_ig"; }

  /// Resume computation after a reduction
  virtual void compute_resume ( Block * block,
				CkReductionMsg * msg) throw(); 

private: // methods

  void compute_resume_ (Block * block, CkReductionMsg * msg) throw();

private: // attributes

  double gamma_;

  double pressure_initial_;

  double solenoidal_fraction_;
  
  double kfmin_;
  double kfmax_;

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

#endif /* ENZO_ENZO_METHOD_MHDTURBULENCEIG_HPP */
