// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoPhysicsFluidProps.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-07-31
/// @brief    [\ref Enzo] Declaration of the EnzoPhysicsGravity class

#ifndef ENZO_GRAVITY_ENZO_PHYSICS_GRAVITY_HPP
#define ENZO_GRAVITY_ENZO_PHYSICS_GRAVITY_HPP

class EnzoPhysicsGravity : public Physics {

  /// @class    EnzoPhysicsGravity
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides a central interface for querying the
  ///    gravitational constant.
  ///
  /// Handling of the gravitational constant depends on simulation's nature:
  ///  - In non-cosmological sims the constant can be set to an arbitrary
  ///    positive value (this helps simplify some test-problems).
  ///  - Users aren't allow to specify the constant's value in cosmological
  ///    simulations. This is because cosmological units are defined such that
  ///    ``4*pi*G*rho_bar`` has the value ``1.0``, where ``rho_bar`` is the
  ///    mean physical matter density of the universe.
  ///
  /// It seems like overkill to dedicate a whole class to storing a single
  /// piece of information like the gravitational constant.
  ///  - However, this information DEFINITELY belongs in a Physics class since
  ///    it needs to be accessed in many (mostly-independent) Method classes
  ///  - More importantly, there isn't a single Method class that can store
  ///    this information. We will soon have 2 Methods that implement distinct
  ///    methods for modelling self-gravity. Plus, it's possible to use a
  ///    Method that implements a background potential independently of either
  ///    of those choices.
  ///  - In the future, it might make sense to dedicate a single Physics object
  ///    to group miscellaneous physical values like the Gravitational constant
  ///
  /// @note
  /// Some imporant notes include:
  ///  - instances of this class can be initialized even if Gravity is not
  ///    considered within a given simulation.
  ///  - the data stored in the class is meant to be considered immutable.
  ///  - enzo::cosmology() should return a consistent value over the lifetime
  ///    of this class (in other words, EnzoPhysicsCosmology must be
  ///    initialized first)

public: // interface

  /// Constructor
  ///
  /// @param grav_const_code_units When positive, this specifies the
  ///     gravitational constant in code units. Otherwise, the gravitational
  ///     constant is assumed to be enzo_constants::grav_const when in cgs
  EnzoPhysicsGravity(double grav_const_code_units)
    : Physics()
  {
    if (grav_const_code_units <= 0.0) {
      grav_constant_codeU_ = -1.0;
    } else {
      ASSERT("EnzoPhysicsGravity::EnzoPhysicsGravity",
             "users aren't allowed to specify the gravitational constant in "
             "a cosmological simulation", enzo::cosmology() == nullptr);
      grav_constant_codeU_ = grav_const_code_units;
    }
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoPhysicsGravity);

  /// CHARM++ migration constructor
  EnzoPhysicsGravity(CkMigrateMessage *m)
    : Physics (m),
      grav_constant_codeU_(-1.0)
  { }

  /// Virtual destructor
  virtual ~EnzoPhysicsGravity()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    Physics::pup(p);
    p | grav_constant_codeU_;
  }

  /// override the virtual method specifying the type
  std::string type() const override { return "gravity"; }

  /// Return the user-configurable gravitational constant in cgs
  double grav_constant_cgs() const noexcept
  {
    if (grav_constant_codeU_ <= 0) { return enzo_constants::grav_constant; }

    Units* units = (Units*)enzo::units();
    double factor = units->density() * (units->time() * units->time());
    return grav_constant_codeU_ / factor;
  }

  /// Return the user-configurable gravitational constant in code units.
  ///
  /// @param units Parameter specifying a reference to a Units instance (or an
  ///     instance of a Units subclass). This is optional in a non-cosmological
  ///     simulation (if a Units-instance is needed in this scenario, it will
  ///     be acquired via enzo::units()). However, this argument is required in
  ///     a cosmological simulation.
  ///
  /// @note
  /// The units argument is required during cosmological simulations since the
  /// the code units are comoving, which means the conversion from cgs to
  /// comoving depends on simulation time. It is caller's responsibility, that
  /// the units argument is configured for the correct time. (We explicitly
  /// avoid doing this within this method to avoid altering the simulation's
  /// global state).
  ///
  /// @note
  /// If the user sets to the real-world value, this should correspond to
  /// 6.67384e-8 cm^3 g^-1 s^-2 when converted to cgs.
  double grav_constant_codeU(const Units& units) const noexcept
  {
    // this first case can only arise in non-cosmological sims
    if (grav_constant_codeU_ > 0) { return grav_constant_codeU_; }

    double factor = units.density() * (units.time() * units.time());
    return enzo_constants::grav_constant * factor;
  }

  double grav_constant_codeU() const noexcept
  {
    ASSERT("EnzoPhysicsGravity::grav_constant_codeU()",
           "requires units arg in a cosmo sim", enzo::cosmology()==nullptr);
    if (grav_constant_codeU_ > 0) { return grav_constant_codeU_; }
    Units* units = (Units*)enzo::units();
    ASSERT("EnzoPhysicsGravity::grav_constant_codeU()",
           "enzo::units() returned a nullptr", units != nullptr);
    return grav_constant_codeU(*units);
  }

protected:
  /// Gravitational constant in code units (if set to the real-world value, it
  /// should correspond to 6.67384e-8 cm^3 g^-1 s^-2 when converted to cgs).
  ///
  /// When this isn't positive, we convert from enzo_constants::grav_constant
  double grav_constant_codeU_;
};

#endif /* ENZO_GRAVITY_ENZO_PHYSICS_GRAVITY_HPP */
