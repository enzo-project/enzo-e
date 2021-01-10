// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMHDVlct.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs June 13 2019
/// @brief    [\ref Enzo] Declaration of the VL + CT (van Leer combined with
///           constrained transport) MHD method
///
/// This class makes relies on the following component classses
///
///    EnzoEquationOfState:       Equation of State for Gas
///    EnzoReconstructor:         Reconstructs primitive variables
///    EnzoRiemann:               Solves the Riemann Problem
///    EnzoConstrainedTransport:  Performs constrained transport
///
/// Some notes on implementation
///    - this Method tracks specific total energy (referred to as total_energy)
///    - We categorize quantities as reconstructable and integrable primitives.
///      Nearly every field overlaps between the two cases. Examples of the
///      primitives for adiabatic, ideal gas (without dual energy formalism):
///        - density         (both integrable and reconstructable)
///        - velocity       (both integrable and reconstructable)
///        - pressure        (just reconstructable)
///        - total_energy    (just integrable)
///        - bfield          (both integrable and reconstructable)
///      When using the dual energy formalism there is also:
///        - internal_energy (just integrable)
///
///    Grouping Objects
///    ----------------
///    Implementation relies upon passing around instances of the Grouping
///    object between different object methods. Originally attempted to pass
///    around vectors of field ids, but that doesn't offer nearly as much
///    flexibility. In particular:
///        - It would make it more difficult to add new fields without altering
///          the interface (e.g. Cosmic Ray energy density and flux density)
///        - Relatedly, by sharing common group names between different
///          grouping objects it is easier to load in relevant fields
///          representing different, but related physical quantities, (e.g.
///          Groupings of primitives and fluxes each always have "density" and
///          "momentum" groups)
///
///    Notes about current groupings:
///        - Currently track 11 different groupings:
///            1. primitive quantities
///            2. interface b-fields (longitudinal B-fields centered at
///               interface between cells)
///            3. reconstructed left primitive fields
///            4. reconstructed right primitive fields
///            5. fluxes in the x-direction
///            6. fluxes in the y-direction
///            7. fluxes in the z-direction
///            8. fields used to accumulate the total change in the conserved
///               versions of all integrable quantities (including flux
///               divergence and source terms)
///            9. face centered fields along x,y,z direction that identify
///               which direction is upwind. (known as weight_group)
///           10. edge centered electric fields
///           11. temp integrable quantities (to store values at the half
///               timestep)
///           12. temp interface b-fields (to store values at the half
///               timestep)
///        - several groups only contain 1 field (e.g. the "density" and
///          "total_energy" groups). In effect, they serve as alias names for
///          the fields
///        - This implementation relies on the fact that the fields in a
///          Grouping are sorted alphabetically (e.g. fields in the "velocity"
///          would always be ordered "velocity_x","velocity_y", "velocity_z")  
///        - All fields within Groupings 1, 8, and 11 are all cell-centered.
///        - All face-centered fields with the exception of (temporary)
///          interface b-fields only include space for values on the interior
///          of the grid. The (temporary) interface bfields all include space
///          for values on the exterior of the grid.
///        - All reconstructed fields are not technically registered as
///          cell-centered fields. They are reused to store reconstructed values
///          along different axes. Consequently, they are formally registered
///          as cell-centered fields (to guarantee that they have enough space).
///        - All reconstructable and integrable primitive quantities have
///          groupings named for them in Groupings 1, 3, 4, 11. We don't
///          construct separate groupings due to the high degree of overlap. To
///          help make the control flow of the calculations more apparent,
///          we create the following aliases for each of the groupings for use
///          in the main for-loop:
///             cur_integrable_prim, out_integrable_prim, reconstructable_prim,
///             reconstructable_priml, reconstructable_primr,
///             integrable_priml, integrable_primr

#ifndef ENZO_ENZO_METHOD_VLCT_HPP
#define ENZO_ENZO_METHOD_VLCT_HPP

class EnzoMethodMHDVlct : public Method {

  /// @class    EnzoMethodMHDVlct
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate VL + CT MHD method

public:
  /// This is defined within the scope of EnzoMethodMHDVlct to avoid polluting
  /// the global namespace.
  ///
  /// @note
  /// This must be public so that it can be passed to helper methods by value
  enum bfield_choice {
    no_bfield = 0,         // pure hydrodynamics
    unsafe_const_uniform,  // an unsafe mode where bfields are assumed to be
                           // const (no interface bfields or CT). This is
                           // provided primarily for debugging)
    constrained_transport  // constrained transport (include interface bfields)
  };

public: // interface

  /// Create a new EnzoMethodMHDVlct object
  EnzoMethodMHDVlct(std::string rsolver,
		    std::string half_recon_name,
		    std::string full_recon_name,
		    double gamma, double theta_limiter,
		    double density_floor,
		    double pressure_floor,
		    std::string mhd_choice,
		    bool dual_energy_formalism,
		    double dual_energy_formalism_eta);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodMHDVlct);

  /// Charm++ PUP::able migration constructor
  EnzoMethodMHDVlct (CkMigrateMessage *m)
    : Method (m),
      eos_(nullptr),
      half_dt_recon_(nullptr),
      full_dt_recon_(nullptr),
      riemann_solver_(nullptr),
      integrable_updater_(nullptr),
      mhd_choice_(bfield_choice::no_bfield),
      integrable_field_list_(),
      reconstructable_field_list_(),
      nested_passive_list_()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Delete EnzoMethodMHDVlct object
  ~EnzoMethodMHDVlct();

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "mhd_vlct"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // methods

  /// returns the bfield_choice enum that matches the input string
  bfield_choice parse_bfield_choice_(std::string choice) const noexcept;

  /// Determines the quantities from (FIELD_TABLE) to be reconstructed and
  /// integrated and that will be integrated.
  ///
  /// @param[in]  eos Pointer to the fluid's EquationOfState object.
  /// @param[in]  mhd_choice Encodes how the integrator will handle B-fields
  /// @param[out] integrable_quantities Reference to a vector that get's filled
  ///     by this function with the integrable quantities (matching names in
  ///     FIELD_TABLE) used by the integrator
  /// @param[out] reconstructable_quantities Reference to a vector that get's
  ///     filled by this function with the names of quantities (matching
  ///     names in FIELD_TABLE) that are used by the integrator for
  ///     reconstruction
  static void determine_quantities_
  (const EnzoEquationOfState *eos, bfield_choice mhd_choice,
   std::vector<std::string> &integrable_quantities,
   std::vector<std::string> &reconstructable_quantities) noexcept;

  /// Checks that the mesh size is sufficiently large to handle the given ghost
  /// depth and confirms that the ghost depth is consistent with the
  /// requirements of the reconstructors
  ///
  /// @param[in] block used to determine the current mesh size and ghost depth
  void check_mesh_and_ghost_size_(Block *block) const;

  /// Converts conservative passive scalars (which are originally densities)
  /// to specific form (basically just divide by density)
  ///
  /// @param block holds data to be processed
  /// @param passive_groups A vector that holds the names of the groups of the
  ///     passive scalars that must be converted
  /// @param conserved_passive_scalars A grouping which holds the group names
  ///     found in passive_groups. The fields within each group hold the
  ///     conserved form (density) of the passive scalars
  /// @param primitive_group A grouping which holds groups named after the
  ///     names in passive_groups. The groups holds the field name where the
  ///     specific form (mass fraction) of the passive scalars will be saved.
  /// @param stale_depth indicates the current stale depth
  void compute_specific_passive_scalars_
  (const std::vector<std::vector<std::string>> passive_lists,
   EFlt3DArray& density, EnzoEFltArrayMap& conserved_passive_scalar_map,
   EnzoEFltArrayMap& specific_passive_scalar_map,
   int stale_depth) const noexcept;

  /// Constructs a map containing the field data for each primitive (except for
  /// the passively advected scalars).
  EnzoEFltArrayMap nonpassive_primitive_map_(Block * block) const throw ();
  
  /// Constructs a map containing the field data for each passively advected
  /// scalar (in conserved form).
  EnzoEFltArrayMap conserved_passive_scalar_map_(Block * block) const throw ();

  /// Computes the fluxes along a given dimension, `dim`, and accumulate the
  /// changes to the integrable quantities in `dUcons_group`
  ///
  /// If using the dual energy formalism, this also computes a part of the
  /// internal energy density source term,
  ///    `dt * pressure * (dvx/dx + dvy/dy + dvz/dz)`
  /// (`vx`, `vy`, `vz` are velocity components & scale factor dependence is
  /// omitted), and adds it to the relevent field in `dUcons_group`. More
  /// specifically, it handles the dimensionally split part of the term
  /// involving the derivative along `dim`. It uses the velocity component
  /// along `dim` at the cell-interfaces (estimated by the Riemann Solver) to
  /// compute the derivatives.
  ///
  /// This function should NOT be modified to directly compute any other source
  /// terms unless they similarly have dependence on dimensional quantites
  /// computed in this function AND can be dimensionally split.
  ///
  /// @param block holds data to be processed
  /// @param dim Dimension along which to compute fluxes. Values of 0, 1, and 2
  ///     correspond to the x, y, and z directions, respectively.
  /// @param reconstructable_group contains the fields holding reconstructable
  ///     quantities that are to be reconstructed. (This includes specific
  ///     passive scalars)
  /// @param cur_bfieldi_group holds the current interface magnetic fields
  /// @param priml_group,primr_group contains the fields that will temporary
  ///     hold the reconstructed reconstructable and integrable quantities
  ///     (and specific passive scalars). The relevant fields should are
  ///     formally defined as cell-centered (to allow for reuse), but during
  ///     calculations, they are treated as face-centered (without having
  ///     values on the exterior faces of the block).
  /// @param pressure_name_l,pressure_name_r are the names of the fields
  ///     that temporarily store the left/right pressure values that are
  ///     computed from the reconstructed primitives. The face-centering is
  ///     expected to match fields contained by priml_group, primr_group
  /// @param flux_group holds field names where the calculated fluxes will be
  ///     stored. The relevant fields should be face-centered along the
  ///     dimension of the calculation (without having values on the exterior
  ///     faces of the block)
  /// @param weight_group holds the temporary weight fields. There is a weight
  ///     field for each spatial direction and the field is face-centered along
  ///     that direction (without including values on the exterior faces of the
  ///     block). The weight field corresponding to dim, is modified to
  ///     indicate the upwind direction (this is included to optionally
  ///     implement the weighting scheme used by Athena++ at a later date)
  /// @param dUcons_group holds field names of fields used to accumulate the
  ///     changes to the conserved forms of the integrable quantities (other
  ///     than Bfield if CT is used) and passive scalars. The changes in the
  ///     integrable quantities and passive scalars due to the fluxes and (any
  ///     source terms) computed by these quantities are all ADDED to these
  ///     fields.
  /// @param interface_velocity_name Expected to be the name of a cell-centered
  ///     field that is used to store temporarily store the velocity along
  ///     `dim` at each cell interface during this function call. It is used to
  ///     compute the internal energy source term (and is required if the dual
  ///     energy formalism is in use). A value of "" indicates that no field
  ///     has been allocated.
  /// @param reconstructor the instance of EnzoReconstructor to use to update
  ///     reconstruct the face centered values
  /// @param stale_depth indicates the current stale depth (before performing
  ///     reconstruction)
  ///
  /// @par Note
  /// It might be worth breaking this into 2 functions (where one of them
  /// handles the reconstruction of the fields and calculation of the flux and
  /// the other additionally handles the accumulation of values in flux_group
  /// and calculates any relevant source terms.
  void compute_flux_
  (int dim, double cur_dt, enzo_float cell_width,
   EnzoEFltArrayMap &reconstructable_map,
   EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
   EFlt3DArray &pressure_l, EFlt3DArray &pressure_r,
   EnzoEFltArrayMap &flux_map, EnzoEFltArrayMap &dUcons_map,
   EFlt3DArray *interface_velocity_arr_ptr, EnzoReconstructor &reconstructor,
   EnzoConstrainedTransport *ct_handler, int stale_depth,
   const std::vector<std::vector<std::string>>& passive_lists) const noexcept;

  /// Allocate temporary fields needed for scratch space and store their names
  /// in the corresponding groupings or return the names. Also allocates all
  /// temporary fields for reconstructable primives that don't already exist in
  /// primitive_group_
  ///
  /// @param block holds data to be processed
  /// @param priml_group,primr_group holds temporary field names where the
  ///     left/right reconstructed face-centered reconstructable and integrable
  ///     primitives are stored. These are registerred as cell-centered fields
  ///     so that they are guaranteed to have enough room to serve as
  ///     face-centered along each dimension. Each grouping will be initialized
  ///     with groups named after the groups contained within primitive_group_.
  /// @param pressure_name_l,pressure_name_r will be set equal to the names of
  ///     the temporary fields that will store the left/right pressure values
  ///     that will be computed from the reconstructed values. Like the fields
  ///     within priml_group and primr_group, these fields will be registered
  ///     as face-centered so that they can be used to hold face-centered
  ///     quantities along each dimension. Note: in the case of an adiabatic,
  ///     ideal gas, pressure is a reconstructable quantity and fields are
  ///     included for it in priml_group and primr_group. In that case, those
  ///     field names are also asigned to pressure_name_l and pressure_name_r.
  /// @param interface_velocity_name This will be set equal to the name of the
  ///     field used to temporarily store the velocity normal to the cell
  ///     interface at the interface. Because it's used to compute the internal
  ///     energy source term, the field is only allocated if the the dual
  ///     energy formalism is in use. Otherwise, it will be set to "".
  /// @param xflux_group,yflux_group,zflux_group The groupings of temporary
  ///     face-centered fields that will store the x, y, and z fluxes. These
  ///     fields don't include the faces on the exterior of the grid. Each
  ///     grouping will be setup with groups named after the names stored in
  ///     integrable_group_names_.
  /// @param dUcons_group Grouping of fields used to accumulate the changes to
  ///     the conserved forms of the integrable quantities and passively
  ///     advected scalars. If CT is used, this grouping won't have space to
  ///     store changes in the magnetic fields (that update is handled
  ///     separately).
  /// @param temp_primitive_group will hold temporary fields identical to all
  ///     the fields held by primitive_group_. Basically these fields hold the
  ///     values computed at the half time-step
  void setup_arrays_
  (Block *block, EnzoEFltArrayMap &primitive_map,
   EnzoEFltArrayMap &temp_primitive_map,
   EnzoEFltArrayMap &conserved_passive_scalar_map,
   EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
   EFlt3DArray &pressure_l, EFlt3DArray &pressure_r,
   EnzoEFltArrayMap &xflux_map, EnzoEFltArrayMap &yflux_map,
   EnzoEFltArrayMap &zflux_map, EnzoEFltArrayMap &dUcons_map) noexcept;

protected: // attributes

  /// Pointer to the equation of state of the fluid
  EnzoEquationOfState *eos_;
  /// Pointer to the reconstructor used to reconstruct the fluid during the
  /// first half time-step (usually nearest-neighbor)
  EnzoReconstructor *half_dt_recon_;
  /// Pointer to the reconstructor used to reconstruct the fluid during the
  /// full time-step
  EnzoReconstructor *full_dt_recon_;
  /// Pointer to the Riemann solver
  EnzoRiemann *riemann_solver_;
  /// Pointer to the integrable quantity updater
  EnzoIntegrableUpdate *integrable_updater_;

  /// Indicates how magnetic fields are handled
  bfield_choice mhd_choice_;

  /// Names of the integrable fields (only includes the field names for
  /// actively advected quantities). These also serve as the keys to the
  /// mappings of arrays used in the calculation
  std::vector<std::string> integrable_field_list_;
  /// Names of the reconstructable primitive fields. These also serve as the
  /// keys to the mappings of arrays used in the calculation
  std::vector<std::string> reconstructable_field_list_;

  /// Lazy initializer of the nested list of fields holding passive scalars
  /// The first sublist holds all field names of quantities that are normally
  /// passively advected. Subsequent lists hold sets of names for scalars whose
  /// total mass fraction must be 1 (like species).
  EnzoNestedPassiveScalarFieldList nested_passive_list_;
};

#endif /* ENZO_ENZO_METHOD_VLCT_HPP */
