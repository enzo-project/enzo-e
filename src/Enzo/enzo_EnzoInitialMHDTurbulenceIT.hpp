// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialMHDTurbulenceIT.hpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:30:57 UTC 2014
/// @date     Wed Aug 24 00:30:57 UTC 2018
/// @brief    [\ref Enzo] Initial conditions for MHD turbulence simulations with 

#ifndef ENZO_ENZO_INITIAL_MHDTURBULENCEIT_HPP
#define ENZO_ENZO_INITIAL_MHDTURBULENCEIT_HPP

class EnzoInitialMHDTurbulenceIT : public Initial {

  /// @class    EnzoInitialMHDTurbulenceIT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for turbulence simulations

public: // interface

  /// CHARM++ constructor
  EnzoInitialMHDTurbulenceIT(double density_initial,
			     double bfieldx_initial,
			     double gamma) throw() 
    : Initial(),
      density_initial_(density_initial),
      bfieldx_initial_(bfieldx_initial),
      gamma_(gamma)

  { }
  
  /// Constructor
  EnzoInitialMHDTurbulenceIT(int cycle, double time,
			     double density_initial,
			     double bfieldx_initial,
			     double gamma) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialMHDTurbulenceIT);

  /// CHARM++ migration constructor
  EnzoInitialMHDTurbulenceIT(CkMigrateMessage *m) 
    : Initial (m), 
      density_initial_(0),
      bfieldx_initial_(0),
      gamma_(1.0)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  (   Block * block,   const Hierarchy * hierarchy   ) throw();


private: // attributes

  double density_initial_;
  double bfieldx_initial_;
  double gamma_;

};

#endif /* ENZO_ENZO_INITIAL_MHDTURBULENCEIT_HPP */

