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

  EnzoMethodBackgroundAcceleration() : zero_acceleration_(false),
                                       mx_(0), my_(0), mz_(0),
                                       gx_(0), gy_(0), gz_(0),
                                       xm_(0), ym_(0), zm_(0),
                                       hx_(0), hy_(0), hz_(0)
                                  {};

  /// Destructor
  virtual ~EnzoMethodBackgroundAcceleration() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodBackgroundAcceleration);

  /// Charm++ PUP::able migration Constructor
  EnzoMethodBackgroundAcceleration (CkMigrateMessage *m)
      : Method(m), zero_acceleration_(false), mx_(0), my_(0), mz_(0),
                   gx_(0), gy_(0), gz_(0), xm_(0), ym_(0), zm_(0),
                   hx_(0), hy_(0), hz_(0)
      { }

  /// CHARM++ Pack / Unpack function
//----------------------------------------------------------------------------

  void pup (PUP::er &p)
  {
    // NOTE: Change this function whenever attributes Change

    TRACEPUP;

    Method::pup(p);

    p | zero_acceleration_;
    p | mx_;
    p | my_;
    p | mz_;
    p | gx_;
    p | gy_;
    p | gz_;
    p | xm_;
    p | ym_;
    p | zm_;
    p | hx_;
    p | hy_;
    p | hz_;
    p | G_four_pi_;

  }

  ///

  ///
  virtual void compute (Block *block) throw();

  virtual std::string name () throw()
  { return "background_acceleration"; }

  virtual double timestep (Block * block) const throw();

  // these don't need to be virtual
  virtual void PointMass( enzo_float * ax,
                          enzo_float * ay,
                          enzo_float * az,
                          Particle * particle,
                          const int rank,
                          const enzo_float cosmo_a,
                          const EnzoConfig * enzo_config,
                          const EnzoUnits * enzo_units,
                          const double dt) throw() ;

  virtual void GalaxyModel( enzo_float * ax,
                          enzo_float * ay,
                          enzo_float * az,
                          Particle * particle,
                          const int rank,
                          const enzo_float cosmo_a,
                          const EnzoConfig * enzo_config,
                          const EnzoUnits * enzo_units,
                          const double dt) throw() ;

protected: // methods

  //double timestep_ (Block * block) const throw();
  void compute_ (Block *block) throw();

protected: // attributes

   /// Convenience. Gravitational constant times 4 pi
   bool zero_acceleration_;
   double G_four_pi_;
   int    mx_, my_, mz_, gx_, gy_, gz_;
   double xm_, ym_, zm_, hx_, hy_, hz_;
};

// make a new class here of acceleration models
// to simplify computation a little bit

#endif /*  ENZO_ENZO_METHOD_BACKGROUND_ACCELERATION_HPP */
