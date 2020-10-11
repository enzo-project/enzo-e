// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrableUpdate.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs June 20 2019
/// @brief    [\ref Enzo] Declaration of Enzo's Integrable Update class. This 
/// is responsible for adding flux and sources terms to integrable quantities.

#include <tuple>

#ifndef ENZO_ENZO_INTEGRABLE_UPDATE_HPP
#define ENZO_ENZO_INTEGRABLE_UPDATE_HPP
class EnzoIntegrableUpdate : public PUP::able
{
  /// @class    EnzoIntegrableUpdate
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the updating of (actively and
  ///           passively) advected integrable quantites

public: // interface

  /// Create a new EnzoIntegrableUpdate instance
  EnzoIntegrableUpdate(std::vector<std::string> integrable_groups,
		       bool skip_B_update,
		       std::vector<std::string> passive_groups) throw();

  /// Virtual destructor
  virtual ~EnzoIntegrableUpdate()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoIntegrableUpdate);

  /// CHARM++ migration constructor for PUP::able
  EnzoIntegrableUpdate (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
    p|first_specific_index_;
    p|density_index_;
    p|passive_groups_;
    p|combined_integrable_groups_;
  }

  /// Iterates through all arrays in `dUcons_map` that are registered as
  /// integrable quantities or are specified with `passive_lists`. All elements
  /// in these arrays are set to `value`.
  ///
  /// @param[in,out] dUcons_map The map holding the arrays that are to be
  ///     modified. These arrays are nominally used to accumulate the changes
  ///     to all integrable and passively advected quantites.
  /// @param[in]     value The value array elements are assigned.
  /// @param[in]     passive_lists A list of lists of keys for passively
  ///     advected scalars. In this method, this is effectively concatenated
  ///     into one list of passive scalars.
  void clear_dUcons_map
  (EnzoEFltArrayMap &dUcons_map, enzo_float value,
   const std::vector<std::vector<std::string>> &passive_lists) const noexcept;

  /// Computes the change in (the conserved form of) the integrable and
  /// passively advected quantites due to the flux divergence along dimension
  /// `dim` over the timestep `dt`. These changes are added to the accumulator
  /// arrays contained by `dUcons_map`.
  ///
  /// @param[in]  dim The dimension along which to compute the flux divergence.
  /// @param[in]  dt The current timestep.
  /// @param[in]  cell_width The cell with along dimension `dim`.
  /// @param[in]  flux_map Map of arrays holding the fluxes computed for
  ///     the current timestep. The values of these fields should be stored
  ///     on the cell faces along the `dim` dimension.
  /// @param[out] dUcons_map Map of arrays where the flux divergence is added
  ///     to. If constrained transport is being used, this will not include
  ///     arrays for the magnetic fields.
  /// @param[in]  stale_depth The stale depth at the time of the function call.
  ///     This should match the stale depth at the time the fluxes were
  ///     computed.
  /// @param[in]  passive_lists A list of lists of keys for passively advected
  ///     scalars. Keys must not be duplicated across more than one list. In
  ///     this method, this is effectively concatenated into one list of
  ///     passive scalar keys.
  void accumulate_flux_component
  (int dim, double dt, double cell_width, EnzoEFltArrayMap &flux_map,
   EnzoEFltArrayMap &dUcons_map, int stale_depth,
   const std::vector<std::vector<std::string>> &passive_lists) const noexcept;

  /// adds flux divergence (and source terms) to the initial integrable
  /// quantities and stores the results in `out_integrable_map`
  ///
  /// @param[in] initial_integrable_map Map of arrays holding the values of
  ///     integrable quantities from the start of the timestep. The fluxes
  ///     will be added to fields held in this grouping. This should also
  ///     contain the passive scalars in specific form (as mass fractions) that
  ///     are to be integrated.
  /// @param[in]  dUcons_map Map of arrays where the net changes to the
  ///     integrable quantities and passively advected quantites are stored.
  ///     If constrained transport is being used, this will not include arrays
  ///     for the magnetic fields.
  /// @param[out] out_integrable_map Map of the fields where the updated
  ///     integrable quantities will be stored (This can be a reference to the
  ///     same Map referenced by initial_integrable_map). The updated passively
  ///     advected scalars will NOT be stored here.
  /// @param[out] out_conserved_passive_scalar Map of arrays where the updated
  ///     passive scalar quantities are stored in conserved form (as densities).
  /// @param[in]  eos Pointer to the fluid's equation of state object. When
  ///     applicable used for placing a density floor.
  /// @param[in] stale_depth The stale depth at the time of the function call
  ///     (the stale_depth must be incremented after this function is called)
  /// @param[in]  passive_lists A list of lists of keys for passively advected
  ///     scalars. The first list holds the keys for quantities that are
  ///     normally passively advected. Subsequent lists group together
  ///     collections pf passively advected scalars whose specific value
  ///     (before and after advection) must sum to 1. Keys must not be
  ///     duplicated across more than one list.
  void update_quantities
  (EnzoEFltArrayMap &initial_integrable_map, EnzoEFltArrayMap &dUcons_map,
   EnzoEFltArrayMap &out_integrable_map,
   EnzoEFltArrayMap &out_conserved_passive_scalar,
   EnzoEquationOfState *eos, int stale_depth,
   const std::vector<std::vector<std::string>> &passive_lists) const;

  /// provides a const vector of all integrable and passively advected scalars
  const std::vector<std::string> combined_integrable_groups() const throw()
  { return combined_integrable_groups_; }

private:

  /// Helper method that updates that takes the initial passively advected
  /// scalars in specific form (as mass fractions) and computes the updated
  /// value in conserved form (as mass densities)
  ///
  /// (This should called before the density is updated)
  void update_passive_scalars_
  (EnzoEFltArrayMap &initial_integrable_map, EnzoEFltArrayMap &dUcons_map,
   EnzoEFltArrayMap &out_conserved_passive_scalar, int stale_depth,
   const std::vector<std::vector<std::string>> &passive_lists) const;

  /// Dynamically allocates and constructs an array of instances of EFlt3DArray
  /// that are loaded from `map` using the ordering of keys in integrable_keys_
  /// @param[in] map Map of arrays holding data related to each integrable
  ///   quantities registered in integrable_keys_
  /// @param[in] stale_depth indicates the current stale_depth for the loaded
  ///   quantities.
  EFlt3DArray* load_integrable_quantities_(EnzoEFltArrayMap &map,
                                           int stale_depth) const;

private: // attributes

  /// Holds key names used to load each integrable quantity component from a
  /// a mapping. Keys for conserved quantities are always listed before the
  /// specfic quantities
  std::vector<std::string> integrable_keys_;
  /// index of the first specific quantity listed in integrable_keys_
  std::size_t first_specific_index_;
  /// index of integrable_keys_ that holds the string "density"
  std::size_t density_index_;

  /// Names of the passively advected groups of fields (e.g. colours)
  std::vector<std::string> passive_groups_;

  /// Name of all integrable quantities and passively advected scalars
  std::vector<std::string> combined_integrable_groups_;

};

#endif /* ENZO_ENZO_INTEGRABLE_UPDATE_HPP */
