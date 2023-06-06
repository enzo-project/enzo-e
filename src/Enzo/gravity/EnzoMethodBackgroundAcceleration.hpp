// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodBackgroundAcceleration.hpp
/// @author   Andrew Emerick (emerick@astro.columbia.edu)
/// @date     2018
/// @brief    [\ref Enzo] Declaration of EnzoMethodBackgroundAcceleration
///
/// Add additional, analytic accelerations

#ifndef ENZO_METHOD_BACKGROUND_ACCELERATION
#define ENZO_METHOD_BACKGROUND_ACCELERATION

class EnzoMethodBackgroundAcceleration : public Method {

  /// @class    enzo_EnzoMethodBackgroundAcceleration
  /// @ingroup  Enzo
  ///
  /// @brief

public: // interface

  /// Create a new EnzoMethodBackgroundAcceleration object
  EnzoMethodBackgroundAcceleration(bool zero_acceleration);

  /// Destructor
  virtual ~EnzoMethodBackgroundAcceleration() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodBackgroundAcceleration);

  /// Charm++ PUP::able migration Constructor
  EnzoMethodBackgroundAcceleration (CkMigrateMessage *m)
      : Method(m), zero_acceleration_(false)
      { }

  /// CHARM++ Pack / Unpack function
//----------------------------------------------------------------------------

  void pup (PUP::er &p)
  {
    // NOTE: Change this function whenever attributes Change

    TRACEPUP;

    Method::pup(p);

    p | zero_acceleration_;
    p | G_four_pi_;

  }

  ///

  ///
  virtual void compute (Block *block) throw();

  virtual std::string name () throw()
  { return "background_acceleration"; }

  virtual double timestep (Block * block) throw();

protected: // methods

  void compute_ (Block *block) throw();

protected: // attributes

   /// Convenience. Gravitational constant times 4 pi
   bool zero_acceleration_;
   double G_four_pi_;
};

// make a new class here of acceleration models
// to simplify computation a little bit

#endif /*  ENZO_ENZO_METHOD_BACKGROUND_ACCELERATION_HPP */
