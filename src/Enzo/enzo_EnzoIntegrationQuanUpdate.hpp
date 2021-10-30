// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrationQuanUpdate.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs June 20 2019
/// @brief    [\ref Enzo] Declaration of Enzo's Integration Quantity Update
/// class. This is responsible for adding flux and sources terms to integration
/// quantities.

#include <tuple>

#ifndef ENZO_ENZO_INTEGRATION_QUAN_UPDATE_HPP
#define ENZO_ENZO_INTEGRATION_QUAN_UPDATE_HPP
class EnzoIntegrationQuanUpdate : public PUP::able
{
  /// @class    EnzoIntegrationQuanUpdate
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the updating of (actively and
  ///           passively) advected integration quantites

public: // interface

  /// Create a new EnzoIntegrationQuanUpdate instance
  EnzoIntegrationQuanUpdate
  (const std::vector<std::string> &integration_quantity_keys,
   const bool skip_B_update) throw();

  /// Virtual destructor
  virtual ~EnzoIntegrationQuanUpdate()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoIntegrationQuanUpdate);

  /// CHARM++ migration constructor for PUP::able
  EnzoIntegrationQuanUpdate (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
    p|integration_keys_;
    p|first_specific_index_;
    p|density_index_;
  }

  /// Iterates through all arrays in `dUcons_map` that are pre-registered
  /// integration quantities or are specified with `passive_list`. All elements
  /// in these arrays are set to `value`.
  ///
  /// @param[in,out] dUcons_map The map holding the arrays that are to be
  ///     modified. These arrays are nominally used to accumulate the changes
  ///     to all integration quantities.
  /// @param[in]     value The value array elements are assigned.
  /// @param[in]     passive_list A list of keys for passive scalars.
  void clear_dUcons_map(EnzoEFltArrayMap &dUcons_map, enzo_float value,
                        const str_vec_t &passive_list) const noexcept;

  /// Computes the change in (the conserved form of) the integration quantities
  /// (including passive scalars) due to the flux divergence along dimension
  /// `dim` over the timestep `dt`. These changes are added to the accumulator
  /// arrays contained by `dUcons_map`.
  ///
  /// @param[in]     dim The dimension along which to compute the flux
  ///     divergence.
  /// @param[in]     dt The current timestep.
  /// @param[in]     cell_width The cell width along dimension `dim`.
  /// @param[in]     flux_map Map of arrays holding the fluxes computed for
  ///     the current timestep. The values of these fields should be stored
  ///     on the cell faces along the `dim` dimension.
  /// @param[in,out] dUcons_map Map of arrays where the flux divergence is added
  ///     to. If constrained transport is being used, this will not include
  ///     arrays for the magnetic fields.
  /// @param[in]     stale_depth The stale depth at the time of the function
  ///     call. This should match the stale depth at the time the fluxes were
  ///     computed.
  /// @param[in]     passive_list A list of keys for passive scalars.
  void accumulate_flux_component(int dim, double dt, enzo_float cell_width,
                                 const EnzoEFltArrayMap &flux_map,
                                 EnzoEFltArrayMap &dUcons_map, int stale_depth,
                                 const str_vec_t &passive_list) const noexcept;

  /// adds flux divergence (and source terms) to the initial integration
  /// quantities and stores the results in `out_integration_map`
  ///
  /// @param[in] initial_integration_map Map of arrays holding the values of
  ///     integration quantities from the start of the timestep. The fluxes
  ///     will be added to fields held in this map. This should also contain
  ///     the passive scalars in conserved form (as densities) that are to be
  ///     integrated.
  /// @param[in]  dUcons_map Map of arrays where the net changes to the
  ///     integration quantities and passively advected quantites are stored.
  ///     If constrained transport is being used, this will not include arrays
  ///     for the magnetic fields.
  /// @param[out] out_integration_map Map of the fields where the updated
  ///     integration quantities will be stored (This can be a reference to the
  ///     same Map referenced by initial_integration_map). The updated
  ///     passively advected scalars will be stored here.
  /// @param[in]  eos Pointer to the fluid's equation of state object. When
  ///     applicable used for placing a density floor.
  /// @param[in]  stale_depth The stale depth at the time of the function call
  ///     (the stale_depth must be incremented after this function is called)
  /// @param[in]  passive_list A list of keys for passive scalars.
  void update_quantities
  (EnzoEFltArrayMap &initial_integration_map,
   const EnzoEFltArrayMap &dUcons_map,
   EnzoEFltArrayMap &out_integration_map, EnzoEquationOfState *eos,
   const int stale_depth, const str_vec_t &passive_list) const;

  /// provides a const vector of all registerred integration keys
  const std::vector<std::string> integration_keys() const throw()
  { return integration_keys_; }

private:

  /// Helper method that updates that takes the initial passively advected
  /// scalars in specific form (as mass fractions) and computes the updated
  /// value in conserved form (as mass densities)
  ///
  /// (This should called before the density is updated)
  void update_passive_scalars_
  (EnzoEFltArrayMap &initial_integration_map,
   const EnzoEFltArrayMap &dUcons_map,
   EnzoEFltArrayMap &out_integration_map, const int stale_depth,
   const str_vec_t &passive_list) const;

  /// Constructs a vector ``EFlt3DArray`` or ``CelloArray<const enzo_float,3>``
  /// that are loaded from `map` using the ordering of keys in integration_keys_
  /// @param[in] map Map of arrays holding data related to each integration
  ///   quantities registered in integration_keys_
  /// @param[in] stale_depth indicates the current stale_depth for the loaded
  ///   quantities.
  const std::vector<EFlt3DArray> load_integration_quan_
  (EnzoEFltArrayMap &map, const int stale_depth) const noexcept;

  const std::vector<CelloArray<const enzo_float, 3>> load_integration_quan_
  (const EnzoEFltArrayMap &map, int stale_depth) const noexcept;

private: // attributes

  /// Holds key names used to load each integration quantity component from a
  /// a mapping. Keys for conserved quantities are always listed before the
  /// specfic quantities. This excludes passively advected scalars.
  std::vector<std::string> integration_keys_;
  /// index of the first specific quantity listed in integration_keys_
  std::size_t first_specific_index_;
  /// index of integration_keys_ that holds the string "density"
  std::size_t density_index_;
};

#endif /* ENZO_ENZO_INTEGRATION_QUAN_UPDATE_HPP */
