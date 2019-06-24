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
///      and uses specific internal energy (refered to as internal_energy)
///    - We categorize quantities as reconstructable and integrable primitives.
///      Nearly every field overlaps between the two cases. Examples of the
///      primitives for an adiabatic gas (without dual energy formalism):
///        - density         (both integrable and reconstructable)
///        - velocitiy       (both integrable and reconstructable)
///        - internal_energy (just reconstructable)
///        - total_energy    (just integrable)
///        - bfield          (both integrable and reconstructable)
///
///    Grouping Objects
///    ----------------
///    Implementation relies upon passing around instances of the Grouping
///    object between different object methods. Originally attempted to pass
///    around vectors of field ids, but that doesn't offer nearly as much
///    flexibility. In particular:
///        - It would make it more difficult to add new fields without altering
///          the interface (e.g. optionally using a dual energy formalism,
///          species, colors, Cosmic Ray energy density and flux density)
///        - Relatedly, by sharing common group names between different
///          grouping objects it is easier to load in relevant fields
///          representing different, but related physical quantities, (e.g.
///          Groupings of primitives and fluxes each always have "density" and
///          "momentum" groups)
///
///    Notes about current groupings:
///        - Currently track 14 different groupings:
///            1. primitive quantities
///            2. interface b-fields (longitudinal B-fields centered at
///               interface between cells)
///            3. reconstructed left primitive fields
///            4. reconstructed right primitive fields
///            5. fluxes in the x-direction
///            6. fluxes in the y-direction
///            7. fluxes in the z-direction
///            8. face centered fields along x,y,z direction that identify
///               which direction is upwind. (known as weight_group)
///            9. edge centered electric fields
///           10. temp integrable quantities (to store values at the half
///               timestep)
///           11. temp interface b-fields (to store values at the half
///               timestep)
///        - several groups only contain 1 field (e.g. the "density" and
///          "total_energy" groups). In effect, they serve as alias names for
///          the fields
///        - This implementation relies on the fact that the fields in a
///          Grouping are sorted alphabetically (e.g. fields in the "velocity"
///          would always be ordered "velocity_x","velocity_y", "velocity_z")  
///        - All fields within the primitive quantites, conserved quanties, and
///          temp conserved quantities are all cell-centered.
///        - All face-centered fields with the exception of (temporary)
///          interface b-fields only include space for values on the interior
///          of the grid. The (temporary) interface bfields all include space
///          for values on the exterior of the grid.
///        - All reconstructed fields are not technically registered as
///          cell-centered fields. They are reused to store reconstructed values
///          along different axes. Consequently, they are formally registerred
///          as cell-centered fields (to guaruntee that they have enough space).
///        - All reconstructable and integrable primitive quantities have
///          groupings named for them in Groupings 1, 3, 4, 10. We don't
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

public: // interface

  /// Create a new EnzoMethodMHDVlct object
  EnzoMethodMHDVlct(std::string rsolver,
		    std::string half_recon_name,
		    std::string full_recon_name,
		    double gamma, double density_floor,
		    double pressure_floor);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodMHDVlct);

  /// Charm++ PUP::able migration constructor
  EnzoMethodMHDVlct (CkMigrateMessage *m)
    : Method (m),
      primitive_group_(NULL),
      bfieldi_group_(NULL),
      eos_(NULL),
      half_dt_recon_(NULL),
      full_dt_recon_(NULL),
      riemann_solver_(NULL),
      integrable_updater_(NULL),
      reconstructable_group_names_(),
      integrable_group_names_(),
      passive_group_names_()
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

  /// Determines the quantities from (FIELD_TABLE) to be reconstructed and
  /// integrated and a list of group names that may or may not include passive
  /// scalars that will be integrated.
  ///
  /// @param eos Pointer to the fluid's EquationOfState object. This is used to
  ///     help determine which quantities are used by the integrator
  /// @param integrable_quantities Reference to a vector that get's filled by
  ///     this function with the integrable quantities (matching names in
  ///     FIELD_TABLE) used by the integrator
  /// @param reconstructable_quantities Reference to a vector that get's filled
  ///     by this function with the reconstructable quantities (matching names
  ///     in FIELD_TABLE) used by the integrator
  /// @param passive_groups Reference to a vector that get's filled by this
  ///     function with the names of groups of passively advected scalars that
  ///     the integrator will advect
  void determine_quantities_
  (EnzoEquationOfState *eos,
   std::vector<std::string> &integrable_quantities,
   std::vector<std::string> &reconstructable_quantities,
   std::vector<std::string> &passive_groups);

  /// Prepare the main groupings used by the integrator
  void setup_groupings_(std::vector<std::string> &integrable_groups,
			std::vector<std::string> &reconstructable_groups,
			std::vector<std::string> &passive_groups);

  /// Checks that the mesh size is big enough given the ghost depth and checks
  /// the ghost depths given the reconstructors
  ///
  /// @param block used to determine the current mesh size and ghost depth
  void check_mesh_and_ghost_size_(Block *block);

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
  ///
  /// @note The names of the fields that hold the conserved form of the passive
  /// scalars are stored in the grouping held by the instance of FieldDescr
  /// that can be instantiated with cello::field_descr()
  void compute_specific_passive_scalars_
  (Block *block, const std::vector<std::string> passive_groups,
   const Grouping conserved_passive_scalars,
   const Grouping primitive_group)
  { }

  /// Computes the fluxes along a given dimension
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
  /// @param reconstructor the instance of EnzoReconstructor to use to update
  ///     reconstruct the face centered values
  /// @param stale_depth indicates the current stale depth (before performing
  ///     reconstruction)
  void compute_flux_(Block *block, int dim,
		     Grouping &reconstructable_group,
		     Grouping &cur_bfieldi_group,
		     Grouping &priml_group, Grouping &primr_group,
		     std::string pressure_name_l,  std::string pressure_name_r,
		     Grouping &flux_group, Grouping &weight_group,
		     EnzoReconstructor &reconstructor, int stale_depth);

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
  ///     quantities along each dimension.
  /// @param xflux_group,yflux_group,zflux_group The groupings of temporary
  ///     face-centered fields that will store the x, y, and z fluxes. These
  ///     fields don't include the faces on the exterior of the grid. Each
  ///     grouping will be setup with groups named after the names stored in
  ///     integrable_group_names_.
  /// @param efield_group Grouping of temporary edge-centered fields where the
  ///     calculated electric fields get stored. Upon initialization, this will
  ///     have one group called "efield". The group will store temporary fields
  ///     to store the electric field along the x, y, and z components. The
  ///     temporary field will be cell centered along the dimension of the
  ///     component and face-centered along the other components, excluding
  ///     exterior faces of the mesh (e.g. The x-component will be
  ///     cell-centered along the x-direction, but face-centered along the y-
  ///     and z- axes)
  /// @param center_efield_name will be set equal to the names of the temporary
  ///     fields where a component of cell-centered E-fields will be
  ///     temporarily stored. This field will get reused for each dimension.
  /// @param weight_group will hold temporary fields inside of a single group
  ///     called "weight". "weight" holds three fields, one for each spatial
  ///     direction. The weight field for a given direction is face-centered
  ///     along that direction (but the field exclude exterior faces of the
  ///     grid) and keeps track of the upwind direction.
  /// @param temp_primitive_group will hold temporary fields identical to all
  ///     the fields held by primitive_group_. Basically these fields hold the
  ///     values computed at the half time-step
  /// @param temp_bfieldi_group will hold temporary fields identical to all
  ///     the fields held by bfieldi_group_. These will hold the values of the
  ///     face-centered fields computed at the half time-step
  void allocate_temp_fields_(Block *block,
			     Grouping &priml_group,
			     Grouping &primr_group,
			     std::string &pressure_name_l,
			     std::string &pressure_name_r,
			     Grouping &xflux_group,
			     Grouping &yflux_group,
			     Grouping &zflux_group,
			     Grouping &efield_group,
			     std::string &center_efield_name,
			     Grouping &weight_group,
			     Grouping &temp_primitive_group,
			     Grouping &temp_bfieldi_group);

  /// Deallocates the temporary fields used for scratch space
  void deallocate_temp_fields_(Block *block, Grouping &priml_group,
			       Grouping &primr_group,
			       std::string pressure_name_l,
			       std::string pressure_name_r,
			       Grouping &xflux_group,
			       Grouping &yflux_group,
			       Grouping &zflux_group,
			       Grouping &efield_group,
			       std::string center_efield_name,
			       Grouping &weight_group,
			       Grouping &temp_primitive_group,
			       Grouping &temp_bfieldi_group);

protected: // attributes

  /// Pointer to a grouping that holds groups named for each of the
  /// reconstructable and integrable quantities. Groups named for scalars hold
  /// 1 field and Groups named for vectors hold 3 quantities (for each
  /// spatial component). All fields related to integrable quantites are
  /// permanent. All fields exclusively related to reconstructable quantities
  /// can be permanent or temporary. Also holds groups named after all known
  /// passive scalar groups. Within a passive scalar group, there are temporary
  /// fields that will be used to hold the passive scalars in specific form.
  Grouping *primitive_group_;

  /// Grouping that holds a group called "bfield". Within that group, there are
  /// three fields. Each field holds a component of the bfield that is
  /// face-centered along the dimension of the component (hence they are
  /// interface bfields)
  Grouping *bfieldi_group_;

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

  /// Names of the reconstructable primitive quantities
  std::vector<std::string> reconstructable_group_names_;
  /// Names of the integrable primitive quantities
  std::vector<std::string> integrable_group_names_;
  /// Names of the groups of passively advected scalars
  std::vector<std::string> passive_group_names_;

};

#endif /* ENZO_ENZO_METHOD_VLCT_HPP */
