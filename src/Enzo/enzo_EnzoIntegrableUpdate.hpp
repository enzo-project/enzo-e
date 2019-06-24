// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoIntegrableUpdate.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs June 20 2019
/// @brief    [\ref Enzo] Declaration of Enzo's Integrable Update class. This 
/// is responsible for adding flux and sources terms to integrable quantities.

// The main reason this is defined is to allow for easy reuse between different
// methods. (This makes use of the EnzoAdvectionFieldLUT struct)

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
		       bool dual_energy_formalism,
		       std::vector<std::string> passive_groups) throw()
  {
    integrable_groups_ = integrable_groups;
    passive_groups_ = passive_groups;
    flagged_groups_ = std::vector<std::string>();

    if (skip_B_update &&
	std::find(integrable_groups_.begin(), integrable_groups_.end(),
		  "bfield") != integrable_groups.end()){
      flagged_groups_.push_back("bfield");
    }

    ASSERT("EnzoIntegrableUpdate",
	   "Not Presently equipped to handle dual_energy_formalism",
	   !dual_energy_formalism);

    setup_lut_();
  }

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
    p|integrable_groups_;
    p|flagged_groups_;
    if (p.isUnpacking()){
      // avoiding PUPing lookup table
      setup_lut_();
    }
    p|passive_groups_;
  }

  /// adds flux divergence (and source terms) to the initial integrable
  /// quantities and stores the results in out_integrable_group
  ///
  /// @param block holds data to be processed
  /// @param initial_integrable_group contains the fields holding the
  ///     the integrable quantities from the start of the timestep. The fluxes
  ///     will be added to fields held in this grouping. This should contain
  ///     the passive scalars in specific form (as mass fractions) that are to
  ///     be integrated.
  /// @param xflux_group,yflux_group,zflux_group contains the fields holding
  ///     the fluxes computed along the x, y, and z directions, respectively.
  /// @param out_integrable_group contains the fields where the updated
  ///     integrable quantities will be stored (This can be passed the same
  ///     object as initial_integrable_group). The updated passively advected
  ///     scalars will NOT be stored here.
  /// @param out_conserved_passive_scalar contains the fields where the updated
  ///     passive scalar quantities are stored in conserved form (as densities).
  /// @param eos Pointer to the fluid's equation of state object (primarily for
  ///     placing density and pressure floor)
  /// @param dt The time time-step overwhich to apply the fluxes
  /// @param stale_depth The stale depth at the time of the function call
  ///     (after the function call, the stale_depth will need to be
  ///     incremented)
  ///
  /// @note Once we start needing to add source terms, it may make sense to
  /// make the consolidation of source terms and fluxes into a separate step
  void update_quantities(Block *block, Grouping &initial_integrable_group,
			 Grouping &xflux_group,
			 Grouping &yflux_group,
			 Grouping &zflux_group,
			 Grouping &out_integrable_group,
			 Grouping &out_conserved_passive_scalar,
			 EnzoEquationOfState *eos, double dt, int stale_depth);

private:

  /// Helper method that updates that takes the initial passively advected
  /// scalars in specific form (as mass fractions) and computes the updated
  /// value in conserved form (as mass densities)
  ///
  /// (This should called before the density is updated)
  void update_passive_scalars_(Block *block, Grouping &initial_integrable_group,
			       Grouping &xflux_group,
			       Grouping &yflux_group,
			       Grouping &zflux_group,
			       Grouping &out_conserved_passive_scalar,
			       double dt, int stale_depth);

private: //methods

  /// Helper function that simply sets up the lookup table
  void setup_lut_();

private: // attributes
  
  /// Names of the quantities to advect
  std::vector<std::string> integrable_groups_;

  /// Names of the quantities that don't need to be updated, but are useful
  /// to be able to access (e.g. Bfields in the case of Constrained Transport)
  std::vector<std::string> flagged_groups_;

  /// struct lookup-table that maps integrable primitive field names to indices
  EnzoAdvectionFieldLUT lut_;

  /// Number of integrable primitive fields in lut_
  int nfields_;

  /// Indices used to iterate over conserved, specific, other categorized
  /// quantities and stop
  int conserved_start_;
  int conserved_stop_;
  int specific_start_;
  int specific_stop_;
  int other_start_;
  int other_stop_;

  /// Names of the passively advected groups of fields (e.g. colours)
  std::vector<std::string> passive_groups_;

};

#endif /* ENZO_ENZO_INTEGRABLE_UPDATE_HPP */
