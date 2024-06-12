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

#include <memory>

// classes for initializing initial conditions
class ScalarInit;
class VectorInit;

enum class InitializerForm{ conserved, primitive };

struct HydroInitPack{

  HydroInitPack() = default;

  HydroInitPack(ScalarInit* density_init, VectorInit* velocity_init,
                ScalarInit* pressure_init, InitializerForm initializer_form);

  /// initializer for density
  std::shared_ptr<ScalarInit> density_init;
  /// initializer for momentum or velocity
  std::shared_ptr<VectorInit> velocity_init;
  /// initializer for pressure or total energy density
  std::shared_ptr<ScalarInit> pressure_init;

  /// specifies the exact quantities that the initializers correspond to:
  ///   - when equal to InitializerForm::primitive, velocity_init and
  ///     pressure_init hold initializers for velocity and thermal pressure
  ///   - when equal to InitializerForm::conserved, velocity_init and
  ///     pressure_init hold initializers for momentum and total energy density
  InitializerForm initializer_form;
};

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
  EnzoInitialInclinedWave(int cycle, double time, ParameterGroup p) throw();

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
      parallel_vel_(std::numeric_limits<double>::min()),
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
    p | parallel_vel_;
    p | pos_vel_;
    p | wave_type_;
  }
public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // functions

  /// prepares initializers for inclined linear jeans wave
  HydroInitPack prepare_jeans_initializers_(bool is_root_block) const noexcept;

  /// prepares initializers for inclined HD waves
  HydroInitPack prepare_HD_initializers_(bool is_root_block) const noexcept;

  /// prepares initializers for inclined MHD waves
  HydroInitPack prepare_MHD_initializers_(VectorInit **a_init) const noexcept;

  /// handles the allocation of initializer of HD quantities for linear waves
  HydroInitPack build_linear_HD_inits_(double density_back, double density_ev,
                                       double etot_back, double etot_ev,
                                       double mom0_back, double mom1_back,
                                       double mom2_back, double mom0_ev,
                                       double mom1_ev, double mom2_ev) const;

  /// Return a vector of known magnetohydrodynamical waves
  std::vector<std::string> mhd_waves_() const throw();

  /// Return a vector of known hydrodynamical waves
  std::vector<std::string> hd_waves_() const throw();

  /// Returns whether parallel_vel_ has been specified
  bool specified_parallel_vel_() const throw()
  { return parallel_vel_ != std::numeric_limits<double>::min(); }

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

  /// denotes a background velocity along the direction that the wave
  /// propagates. A value of std::numeric_limits<double>::min() means that
  /// it hasn't been set (in which case default values are employed). This can
  /// only be set for hydrodynamical waves
  double parallel_vel_;

  /// whether the wave speed is in the positive direction. For purely
  /// hydrodynamic entropy waves, this has no effect if parallel_vel_ is set.
  bool pos_vel_;

  /// Determines initial values of the wave
  std::string wave_type_;
};

#endif /* ENZO_ENZO_INITIAL_INCLINED_WAVE_HPP */
