// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoPhysicsFluidProps.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-04-05
/// @brief    [\ref Enzo] Declaration of the EnzoPhysicsFluidProps class

#ifndef ENZO_ENZO_PHYSICS_FLUID_PROPS_HPP
#define ENZO_ENZO_PHYSICS_FLUID_PROPS_HPP

class EnzoPhysicsFluidProps : public Physics {

  /// @class    EnzoPhysicsFluidProps
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides a central interface for querying the
  ///    physical properties of fluids.
  ///
  /// The data stored in the class is meant to be considered immutable.

public: // interface

  /// Constructor - this needs to be refactored in the future.
  EnzoPhysicsFluidProps() noexcept;

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoPhysicsFluidProps);

  /// CHARM++ migration constructor
  EnzoPhysicsFluidProps(CkMigrateMessage *m)
    : Physics (m)
  { }

  /// Virtual destructor
  virtual ~EnzoPhysicsFluidProps()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p|de_config_;
    p|fluid_floor_config_;
    p|gamma_;
    p|mol_weight_;
  }

  /// Access the dual energy configuration
  const EnzoDualEnergyConfig& dual_energy_config() const noexcept
  { return de_config_; }

  const EnzoFluidFloorConfig& fluid_floor_config() const noexcept
  { return fluid_floor_config_; }

  enzo_float gamma() const noexcept
  { return gamma_; }

  enzo_float mol_weight() const noexcept
  { return mol_weight_; }

public: // virtual methods
  virtual std::string type() const { return "fluid_props"; }

private: // attributes

  // NOTE: change pup() function whenever attributes change
  EnzoDualEnergyConfig de_config_;
  EnzoFluidFloorConfig fluid_floor_config_;

  /// the adiabatic index used for hydro-related calculations.
  ///
  /// When Grackle is configured with primordial_chemistry>1, its routines will
  /// use a different value of gamma.
  enzo_float gamma_;

  /// the nominal mean molecular weight (i.e. this multiplied by the mass of
  /// hydrogen gives mean molecular mass)
  ///
  /// this is largely ignored when Grackle is in use
  enzo_float mol_weight_;
};

#endif /* ENZO_ENZO_PHYSICS_FLUID_PROPS_HPP */
