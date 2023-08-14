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

  /// @class    EnzoMethodBackgroundAcceleration
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Evaluates an analytic background potential

public: // interface

  /// Create a new EnzoMethodBackgroundAcceleration object
  EnzoMethodBackgroundAcceleration(ParameterAccessor &p);

  /// Destructor
  virtual ~EnzoMethodBackgroundAcceleration() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodBackgroundAcceleration);

  /// Charm++ PUP::able migration Constructor
  EnzoMethodBackgroundAcceleration (CkMigrateMessage *m)
    : Method(m),
      zero_acceleration_(false),
      G_four_pi_(0.0),
      potential_center_xyz_{}, // fills array with zeros
      flavor_(""),
      galaxy_pack_dfltU_(nullptr),
      point_mass_pack_dfltU_(nullptr)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: Change this function whenever attributes Change

    TRACEPUP;

    Method::pup(p);

    p | zero_acceleration_;
    p | G_four_pi_;
    p | potential_center_xyz_;
    p | flavor_;
    p | galaxy_pack_dfltU_;
    p | point_mass_pack_dfltU_;
  }

  ///
  virtual void compute (Block *block) throw();

  virtual std::string name () throw()
  { return "background_acceleration"; }

  virtual double timestep (Block * block) throw();

  const EnzoPotentialConfigGalaxy* try_get_config_galaxy() const
  { return galaxy_pack_dfltU_.get(); }

protected: // methods

  void compute_ (Block *block) throw();

protected: // attributes

  /// Convenience. Gravitational constant times 4 pi
  bool zero_acceleration_;
  double G_four_pi_;

  /// location of the center of the potential, in code units
  std::array<double,3> potential_center_xyz_;

  /// the type of background potential to be used
  std::string flavor_;

  /// stores the parameter-pack for a GalaxyPotential in default units
  ///
  /// (either this or point_mass_pack_dfltU_ must be non-null -- but not both)
  std::unique_ptr<EnzoPotentialConfigGalaxy> galaxy_pack_dfltU_;

  /// stores the parameter-pack for a PointMass in default units
  ///
  /// (either this or galaxy_pack_dfltU_ must be non-null -- but not both)
  std::unique_ptr<EnzoPotentialConfigPointMass> point_mass_pack_dfltU_;
};


#endif /*  ENZO_ENZO_METHOD_BACKGROUND_ACCELERATION_HPP */
