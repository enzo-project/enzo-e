#ifndef ENZO_ENZO_INITIAL_LINEAR_WAVE_HPP
#define ENZO_ENZO_INITIAL_LINEAR_WAVE_HPP

class EnzoInitialLinearWave : public Initial {
  /// @class    EnzoInitialLinearWave
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the Linear Wave test problem for
  /// the VLCT method

public: // interface
  
  /// Constructor
  EnzoInitialLinearWave(int cycle, double time, double alpha, double beta,
			double gamma, int wave_type) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialLinearWave);

  /// CHARM++ migration constructor
  EnzoInitialLinearWave(CkMigrateMessage *m) 
    : Initial (m),
      alpha_(0.0),
      beta_(0.0),
      gamma_(0.0),
      wave_type_(0)
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


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// alpha and beta are angles defined between the axis along which the 1D
  /// wave is defined, and the x,y,z coordinate system
  double alpha_;
  double beta_;

  /// Ideal gas law constant
  double gamma_;

  /// Determines initial values of the wave
  int wave_type_;
};

#endif /* ENZO_ENZO_INITIAL_Linear_Wave_HPP */
