// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialBurkertBodenheimer.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-04-03
/// @brief    [\ref Enzo] 3D array of Collapsing spheres

#ifndef ENZO_ENZO_INITIAL_BURKERTBODENHEIMER_HPP
#define ENZO_ENZO_INITIAL_BURKERTBODENHEIMER_HPP

class EnzoInitialBurkertBodenheimer : public Initial {

  /// @class    EnzoInitialBurkertBodenheimer
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] initial conditions for 3D array of collapsing spheres

public: // interface

  /// CHARM++ constructor
  EnzoInitialBurkertBodenheimer
  (int cycle, double time,
   int rank,
   const int array[3],
   double radius_relative,
   double particle_ratio,
   double mass,
   double temperature,
   int densityprofile
   ) throw()
    : Initial (cycle,time),
      rank_(rank),
      radius_relative_(radius_relative),
      particle_ratio_(particle_ratio),
      mass_(mass),
      temperature_(temperature),
      densityprofile_(densityprofile)
  {
    array_[0] = array[0];
    array_[1] = array[1];
    array_[2] = array[2];
  }    
  
  /// Constructor
  EnzoInitialBurkertBodenheimer(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialBurkertBodenheimer);

  /// CHARM++ migration constructor
  EnzoInitialBurkertBodenheimer(CkMigrateMessage *m)
    : Initial (m),
      rank_(0),
      radius_relative_(0.0),
      particle_ratio_(0.0),
      mass_(0.0),
      temperature_(0.0),
      densityprofile_(1)
  {
    array_[0] = 0;
    array_[1] = 0; 
    array_[2] = 0;
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // attributes

  /// Rank of the problem
  int rank_;
  
  /// Size of the array of Collapsing spheres
  int array_[3];

  /// Relative radius
  double radius_relative_;

  /// Relative ratio of mass distribution between particles and fields
  // (0.0 no particles, 1.0 no fields)
  double particle_ratio_;

  /// Mass of a sphere (ellipsoid)
  double mass_;

  /// Temperature
  double temperature_;

  /// Density Profile
  /// 1 = Uniform Density
  /// 2 = 1/r^2 power law profile 
  int densityprofile_;
};

#endif /* ENZO_ENZO_INITIAL_BURKERTBODENHEIMER_HPP */

