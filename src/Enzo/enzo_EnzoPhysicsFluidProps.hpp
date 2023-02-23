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

  /// Constructor
  EnzoPhysicsFluidProps(const EnzoDualEnergyConfig& de_config,
                        const EnzoFluidFloorConfig& fluid_floor_config,
                        enzo_float gamma, enzo_float mol_weight) noexcept
  : Physics(),
    de_config_(de_config),
    fluid_floor_config_(fluid_floor_config),
    gamma_(gamma),
    mol_weight_(mol_weight)
  { }

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

  // the following methods are mostly just methods of this class because it
  // isn't obvious where else to put them

  /// applies the pressure floor to the specific total energy field. If using
  /// the dual-energy formalism, it is also applied to the internal energy
  /// and it synchronize the internal energy and total energy fields. If the
  /// EOS is barotropic, this does nothing.
  ///
  /// @param[in,out] integration_map Map holding integration quantities that
  ///     will be used to apply the floor. It must also include a
  ///     "total_energy" entry (unless the EOS is barotropic) upon which the
  ///     floor is applied.
  /// @param[in]     stale_depth indicates the current stale_depth for the
  ///     supplied cell-centered quantities
  ///
  /// Unlike the initial conception of the dual-energy formalism (or the
  /// version used in Enzo's ppm integrator) this assumes that synchronization
  /// is a local operation that doesn't require data about neighboring cells
  /// (similar to the implementation of the dual energy formalsim in Enzo's
  /// Runge Kutta and MHD with Constrained Transport solvers).
  ///
  /// @note
  /// It's unclear if this really belongs as a method of EnzoPhysicsFluidProps
  void apply_floor_to_energy_and_sync(EnzoEFltArrayMap &integration_map,
                                      const int stale_depth) const;

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
