// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrableUpdate.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs June 20 2019
/// @brief    [\ref Enzo] Declaration of Enzo's Integrable Update class. This 
/// is responsible for adding flux and sources terms to integrable quantities.

// The main reason this is defined is to allow for easy reuse between different
// methods. (This makes use of the EnzoAdvectionFieldLUT struct)

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
    p|integrable_grouping_items_;
    p|first_specific_index_;
    p|density_index_;
    p|passive_groups_;
    p|combined_integrable_groups_;
  }

  /// Iterates through all fields in `dUcons_group` that are included in groups
  /// registered as integrable quantities or passive scalars (but are not
  /// flagged) and set all of their elements set equal to `value`.
  ///
  /// @param block holds data to be processed
  /// @param dUcons_group The grouping of fields to be modified. These fields
  ///     are nominally used to accumulate the changes to all integrable and
  ///     passively advected quantites.
  /// @param value the value to assign to all elements of the fields
  void clear_dUcons_group(Block *block, Grouping &dUcons_group,
			  enzo_float value) const;

  /// Computes the change in (the conserved form of) the integrable and
  /// passively advected quantites due to the flux divergence along dimension
  /// `dim` over the timestep `dt`. These changes are added to the accumulator
  /// fields specified by `dUcons_group`.
  ///
  /// @param block holds data to be processed
  /// @param dim The dimension along which to compute the flux divergence
  /// @param dt The current timestep.
  /// @param flux_group Grouping of fields that hold the fluxes computed for
  ///     the current timestep. The values of these fields should be stored
  ///     on the cell faces along the `dim` dimension.
  /// @param dUcons_group The grouping of fields that the flux divergence is
  ///     added to.
  /// @param stale_depth The stale depth at the time of the function call.
  ///     This should match the stale depth at the time the fluxes were
  ///     computed.
  void accumulate_flux_component(Block *block, int dim, double dt,
				 Grouping &flux_group, Grouping &dUcons_group,
				 int stale_depth) const;

  /// adds flux divergence (and source terms) to the initial integrable
  /// quantities and stores the results in out_integrable_group
  ///
  /// @param block holds data to be processed
  /// @param initial_integrable_group contains the fields holding the
  ///     the integrable quantities from the start of the timestep. The fluxes
  ///     will be added to fields held in this grouping. This should contain
  ///     the passive scalars in specific form (as mass fractions) that are to
  ///     be integrated.
  /// @param dUcons_group contains the fields where the net changes to the
  ///     integrable quantities and passively advected quantites are stored.
  ///     If CT is being used, this will not include changes to the magnetic
  ///     fields.
  /// @param out_integrable_group contains the fields where the updated
  ///     integrable quantities will be stored (This can be passed the same
  ///     object as initial_integrable_group). The updated passively advected
  ///     scalars will NOT be stored here.
  /// @param out_conserved_passive_scalar contains the fields where the updated
  ///     passive scalar quantities are stored in conserved form (as densities).
  /// @param eos Pointer to the fluid's equation of state object. When
  ///     applicable used for placing density and pressure floors and
  ///     synchronizing the internal energy with the total energy
  /// @param stale_depth The stale depth at the time of the function call
  ///     (the stale_depth needs to be incremented after this function is
  ///     called)
  void update_quantities(Block *block, Grouping &initial_integrable_group,
			 Grouping &dUcons_group,
			 Grouping &out_integrable_group,
			 Grouping &out_conserved_passive_scalar,
			 EnzoEquationOfState *eos, int stale_depth) const;

  /// provides a const vector of all integrable and passively advected scalars
  const std::vector<std::string> combined_integrable_groups() const throw()
  { return combined_integrable_groups_; }

private:

  /// Helper method that updates that takes the initial passively advected
  /// scalars in specific form (as mass fractions) and computes the updated
  /// value in conserved form (as mass densities)
  ///
  /// (This should called before the density is updated)
  void update_passive_scalars_(Block *block, Grouping &initial_integrable_group,
			       Grouping &dUcons_group,
			       Grouping &out_conserved_passive_scalar,
			       int stale_depth) const;

  /// Constructs an array of instances of EFlt3DArray based on the data stored
  /// within combined_integrable_groups_
  /// @param block The mesh block from which data will be loaded
  /// @param grouping Reference to a grouping that includes groups named for
  ///   all registerred quantities. The data is actually loaded from the fields
  ///   in block with names matching the relevant field names specified in this
  ///   object.
  /// @param stale_depth indicates the current stale_depth for the loaded
  ///   quanties.
  EFlt3DArray* load_integrable_quantites_(Block *block, Grouping & grouping,
					  int stale_depth) const;

private: // attributes

  /// Holds group-index pairs used to load each field representing
  /// integrable quantities from a Grouping object. Conserved quantities are
  /// always listed before the specific quantites.
  std::vector<std::pair<std::string,int>> integrable_grouping_items_;
  /// first index of integrable_grouping_items_ corresponding to specific
  /// quantities 
  std::size_t first_specific_index_;
  /// index of integrable_grouping_items_ corresponding to density
  std::size_t density_index_;

  /// Names of the passively advected groups of fields (e.g. colours)
  std::vector<std::string> passive_groups_;

  /// Name of all integrable quantities and passively advected scalars
  std::vector<std::string> combined_integrable_groups_;

};

#endif /* ENZO_ENZO_INTEGRABLE_UPDATE_HPP */
