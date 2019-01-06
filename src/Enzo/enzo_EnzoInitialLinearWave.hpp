#ifndef ENZO_ENZO_INITIAL_LINEAR_WAVE_HPP
#define ENZO_ENZO_INITIAL_LINEAR_WAVE_HPP

// classes for initializing initial conditions
class ScalarInit;
class VectorInit;

// Rotation is responsible for rotating axes.
// This class will be slightly extended and reused for implementing cosmic rays
class Rotation;

class EnzoInitialLinearWave : public Initial {
  /// @class    EnzoInitialLinearWave
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the Linear Wave test problem for
  /// the VLCT method. Assumes adiabatic gas.

public: // interface

  /// Constructor
  EnzoInitialLinearWave(int cycle, double time, double alpha, double beta,
			double gamma, std::string wave_type) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialLinearWave);

  /// CHARM++ migration constructor
  EnzoInitialLinearWave(CkMigrateMessage *m) 
    : Initial (m),
      alpha_(0.0),
      beta_(0.0),
      gamma_(0.0),
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
    p | wave_type_;
  }
public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // functions

  void prepare_initializers_(ScalarInit *density_init,
			     ScalarInit *total_energy_init, 
			     VectorInit *momentum_init,
			     VectorInit *a_init,
			     Rotation &rot);

  bool valid_wave_type_();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// alpha and beta are angles defined between the axis along which the 1D
  /// wave is defined, and the x,y,z coordinate system
  double alpha_;
  double beta_;

  /// Ideal gas law constant
  double gamma_;

  /// Determines initial values of the wave
  std::string wave_type_;
};

#endif /* ENZO_ENZO_INITIAL_Linear_Wave_HPP */
