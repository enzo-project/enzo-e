// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialInclinedWave.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thur April 19 2019
/// @brief    [\ref Enzo] Initialization routine for inclined linear MHD waves
///           and inclined circularly polarized alfven waves detailed in
///           Gardiner & Stone (2008). These are used to test the VL+CT MHD
///           integrator.

#ifndef ENZO_ENZO_INITIAL_INCLINED_WAVE_HPP
#define ENZO_ENZO_INITIAL_INCLINED_WAVE_HPP

// classes for initializing initial conditions
class ScalarInit;
class VectorInit;

// class for rotating axes. Will be (slightly) extended and reused for
// implementing cosmic rays
class Rotation;

class EnzoInitialInclinedWave : public Initial {
  /// @class    EnzoInitialInclinedWave
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the inclined Linear Wave and
  /// Circularly Polarized Alfven wave test problems for the VLCT method.
  /// Assumes adiabatic gas.

public: // interface

  /// Constructor
  EnzoInitialInclinedWave(int cycle, double time, double alpha, double beta,
			double gamma, double amplitude, double lambda,
			bool pos_vel, std::string wave_type) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialInclinedWave);

  /// CHARM++ migration constructor
  EnzoInitialInclinedWave(CkMigrateMessage *m) 
    : Initial (m),
      alpha_(0.0),
      beta_(0.0),
      gamma_(0.0),
      amplitude_(0.0),
      lambda_(0.0),
      pos_vel_(true),
      wave_type_("")
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: update whenever attributes change

    TRACEPUP;

    Initial::pup(p);
    p | alpha_;
    p | beta_;
    p | gamma_;
    p | amplitude_;
    p | lambda_;
    p | pos_vel_;
    p | wave_type_;
  }
public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // functions

  void prepare_initializers_(ScalarInit **density_init,
			     ScalarInit **etot_dens_init, 
			     VectorInit **momentum_init,
			     VectorInit **a_init);

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// alpha and beta are angles defined between the axis along which the 1D
  /// wave is defined, and the x,y,z coordinate system
  double alpha_;
  double beta_;

  /// adiabatic index
  double gamma_;

  /// Amplitude of the perturbation
  double amplitude_;

  /// wavelength
  double lambda_;

  /// If the wave speed is in the positive direction
  bool pos_vel_;

  /// Determines initial values of the wave
  std::string wave_type_;
};

#endif /* ENZO_ENZO_INITIAL_INCLINED_WAVE_HPP */
