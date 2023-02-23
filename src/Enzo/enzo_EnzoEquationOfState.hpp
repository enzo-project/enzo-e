// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEquationOfState.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Equation Of State
/// abstract base class. It will be subclassed to implement specific types of
/// equations of state.

//#define RAISE_FLOOR_ERROR

#ifndef ENZO_ENZO_EQUATIONOFSTATE_HPP
#define ENZO_ENZO_EQUATIONOFSTATE_HPP

// Among it's EOS-related responsbilities, EnzoEquationOfState, is responsible
// for the application of the Dual Energy Formalism, (when specified for
// non-barotropic equations of state). Currently, implementations of the Dual
// Energy Formalism are expected to more closely resemble implementations in
// Enzo's Runge-Kutta and MHD with Constrained Transport integrators. These
// exhibit 3 main differences from the original conception (implemented in
// Enzo's ppm integrator):
//     1. internal energy is always used to compute pressure. In the original
//        conception, pressure could be computed from total energy or
//        internal energy (the decision was independent of synchronization).
//     2. Unlike the original conception, both pressure and internal energy are
//        not reconstructed separately. Implementations are currently expected
//        to just reconstruct pressure and compute internal energy from the
//        reconstructed quantities.
//     3. Synchronization of the total and internal energies is a local
//        operation that doesn't require knowledge cell neighbors. In the
//        original conception, knowledge of the immediate neighbors had been
//        required (thus, each synchronization incremented the stale depth).
//
// To allow synchronization to data from neigboring cells, an additional method
// would be required that indicates the staling_rate of the synchronization.

class EnzoEquationOfState : public PUP::able
{
  /// @class    EnzoEquationOfState
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates equation of state of fluid

public: // interface

  /// Create a new EnzoEquationOfState
  EnzoEquationOfState() throw()
  {}

  ~EnzoEquationOfState()
  { }

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoEquationOfState);

  /// CHARM++ migration constructor for PUP::able
  EnzoEquationOfState (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
  }

  /// returns whether the equation of state is barotropic
  virtual bool is_barotropic() const = 0;

  /// returns adiabatic index - only needs to be a reasonable number of non-
  /// barotropic
  virtual enzo_float get_gamma() const = 0;

  /// returns isothermal sound speed - only needs to be reasonable for a
  /// barotropic EOS
  virtual enzo_float get_isothermal_sound_speed() const = 0;

};

#endif /* ENZO_ENZO_EQUATIONOFSTATE_HPP */
