// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMHDVlct.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs June 13 2019
/// @brief    [\ref Enzo] Declaration of the VL + CT (van Leer combined with
///           constrained transport) MHD method
///
/// This class relies on the following component classses
///
///    EnzoEquationOfState:       Equation of State for Gas
///    EnzoReconstructor:         Reconstructs primitive variables
///    EnzoRiemann:               Solves the Riemann Problem
///    EnzoIntegrableUpdate:      handles the updating of advected quantites
///
/// This class can be run with and without magnetic fields. When run with
/// magnetic fields, it makes use of the following component:
///    EnzoBfieldMethod:  Performs constrained transport
///
/// Some notes on implementation
///    - this Method tracks specific total energy (referred to as total_energy)
///    - We categorize quantities as reconstructable and integrable primitives.
///      Nearly every field overlaps between the two cases. Examples of the
///      primitives for adiabatic, ideal gas (without dual energy formalism):
///        - density         (both integrable and reconstructable)
///        - velocity        (both integrable and reconstructable)
///        - pressure        (just reconstructable)
///        - total_energy    (just integrable)
///      When supporting magnetic fields, there is also:
///        - bfield          (both integrable and reconstructable)
///      When using the dual energy formalism there is also:
///        - internal_energy (just integrable)
///
///    EnzoEFltArrayMap Objects
///    ------------------------
///    The implementation relies upon passing around instances of
///    EnzoEFltArrayMap between the methods of the various components. This
///    class implements a map/dictionary that holds instances of EFlt3DArray.
///    All arrays in a given map are assumed to have the same shape.
///
///    Currently, 8 different maps are used:
///        1. primitive_map: holds the integrable and reconstructed quantities
///           that are stored in Cello fields. It also holds the values of all
///           passively advected scalars, in specific form, which are stored
///           in temporary arrays
///        2. temp_primive_map: This holds temporary arrays for each of the
///           quantities in primitive_map. These arrays are used to store the
///           estimated state at the partial timestep.
///        3. priml_map: holds left reconstructed primitive fields (has the
///           same keys as primitive_map)
///        4. primr_map: holds right reconstructed primitive fields (has the
///           same keys as primitive_map)
///        5. xflux_map: holds fluxes in the x-direction
///        6. yflux_map: holds fluxes in the y-direction
///        7. zflux_map: holds fluxes in the z-direction
///        8. dUcons_map: holds arrays which are used to accumulate the total
///           change in the conserved versions of all integrable quantities
///           (includes both flux divergence and source terms)
///    Note:
///        - All arrays in maps 1, 2, and 8 have the same shapes as
///          cell-centered Cello Fields
///        - All arrays in priml_map and primr_map technically have the shape
///          of a cell-centered field, but they are treated as though they have
///          the shape of a face-centered field. This is done so that they can
///          be reused for different axes.
///        - For the purposes of these enumerated maps, we assume that the
///          length of a face-centered array along the dimension with
///          face-centering is 1 less than that of a cell-centered array
///        - All reconstructed fields are not technically registered as
///          cell-centered fields.  Consequently, they are formally registered
///          as cell-centered fields (to guarantee that they have enough space).
///        - All reconstructable and integrable primitive quantities have
///          key-array pairs named for them in maps number 1, 2, 3, and 4.

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
      bfield_method_(nullptr),
      integrable_field_list_(),
      reconstructable_field_list_(),
      lazy_passive_list_()
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
  void check_mesh_and_ghost_size_(Block *block) const noexcept;

  /// Converts conservative passive scalars (which are originally densities)
  /// to specific form (basically just divide by density)
  ///
  /// @param[in]  passive_list A list of keys for passive scalars.
  /// @param[in]  density Array holding the current density values
  /// @param[in]  conserved_passive_scalar_map Map of the arrays containing the
  ///     current values of the passively advected scalars (in conserved form)
  /// @param[out] specific_passive_scalar_map Map of arrays where the specific
  ///     form of the scalars will be stored.
  /// @param[in]  stale_depth The current stale depth
  void compute_specific_passive_scalars_
  (const str_vec_t &passive_list, EFlt3DArray& density,
   EnzoEFltArrayMap& conserved_passive_scalar_map,
   EnzoEFltArrayMap& specific_passive_scalar_map,
   int stale_depth) const noexcept;

  /// Constructs a map containing the field data for each primitive (except for
  /// the passively advected scalars).
  EnzoEFltArrayMap nonpassive_primitive_map_(Block * block) const noexcept;
  
  /// Constructs a map containing the field data for each passively advected
  /// scalar (in conserved form).
  EnzoEFltArrayMap conserved_passive_scalar_map_(Block * block) const noexcept;

  /// Computes the fluxes along a given dimension, `dim`, and accumulate the
  /// changes to the integrable quantities in `dUcons_map`
  ///
  /// If using the dual energy formalism, this also computes a part of the
  /// internal energy density source term,
  ///    `dt * pressure * (dvx/dx + dvy/dy + dvz/dz)`
  /// (`vx`, `vy`, `vz` are velocity components & scale factor dependence is
  /// omitted), and adds it to the 'internal_energy' entry in `dUcons_map`.
  /// More specifically, it handles the dimensionally split part of the term
  /// involving the derivative along `dim`. The velocity component along `dim`
  /// at the cell-interfaces (estimated by the Riemann Solver) to compute the
  /// derivatives.
  ///
  /// This function should NOT be modified to directly compute any other source
  /// terms unless they similarly have dependence on dimensional quantites
  /// computed in this function AND can be dimensionally split.
  ///
  /// @param[in]     dim Dimension along which to compute fluxes. Values of 0,
  ///     1, and 2 correspond to the x, y, and z directions, respectively.
  /// @param[in]     dt The current timestep.
  /// @param[in]     cell_width The cell width along dimension `dim`.
  /// @param[in]     reconstructable_map Map of arrays holding cell-centered
  ///     primitive quantities that are to be reconstructed (This includes
  ///     specific passive scalars).
  /// @param[in]     priml_map,primr_map Maps of arrays used to temporarily
  ///     hold the left/right reconstructed face-centered reconstructable and
  ///     integrable quantities. These arrays should have the shape of a
  ///     cell-centered field, but are treated as though they have the shape of
  ///     a that is face-centered along `dim` (If a cell-centered field holds
  ///     `N` elements along `dim`, then such a face-centered field should only
  ///     have `N-1` elements along `dim`).
  /// @param[in]     pressure_l,pressure_r Arrays used to temporarily store the
  ///     left/right pressure values. The shape of these arrays should be the
  ///     same as in priml_map,primr_map. Note, pressure_l (pressure_r) is
  ///     allowed to be a shallow copy of an array in priml_map (primr_map).
  /// @param[in]     flux_map Holds arrays where the calculated fluxes
  ///     will be stored. The arrays should be face-centered along `dim`.
  ///     If a cell-centered field holds `N` elements along `dim`, then this
  ///     should only hold `N-1` elements along `dim`.
  /// @param[in,out] dUcons_map Map of arrays where the changes to the
  ///     integrable quantities are accumulated. If constrained transport is
  ///     being used, this won't include arrays for the magnetic fields.
  /// @param[in]     interface_velocity_arr_ptr Pointer to an array to
  ///     temporarily hold the computed component of the velocity at the cell
  ///     interfaces along `dim`. If a cell-centered field holds `N` elements
  ///     along `dim`, then this is only used to store `N-1` elements. This
  ///     quantity is used to compute the internal energy source term (needed
  ///     under the dual energy formalism). If the value is `nullptr`, then the
  ///     interface velocity is not stored in the array.
  /// @param[in]     reconstructor the instance of EnzoReconstructor to use to
  ///     update reconstruct the face centered values
  /// @param[in,out] bfield_method When using running with bfield handling, this
  ///     is a pointer to an instance of EnzoBfieldMethod. During the function
  ///     call, the internal state is updated. If not handling bfields, this
  ///     should be a `nullptr`.
  /// @param[in]     stale_depth indicates the current stale depth (before
  ///     performing reconstruction)
  /// @param[in]     passive_list A list of keys for passively advected scalars.
  ///
  /// @par Note
  /// It might be worth breaking this into 2 functions (where one of them
  /// handles the reconstruction of the fields and calculation of the flux and
  /// the other additionally handles the accumulation of values in flux_map
  /// and calculates any relevant source terms.
  void compute_flux_
  (int dim, double cur_dt, enzo_float cell_width,
   EnzoEFltArrayMap &reconstructable_map,
   EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
   EFlt3DArray &pressure_l, EFlt3DArray &pressure_r,
   EnzoEFltArrayMap &flux_map, EnzoEFltArrayMap &dUcons_map,
   EFlt3DArray *interface_velocity_arr_ptr, EnzoReconstructor &reconstructor,
   EnzoBfieldMethod *bfield_method, int stale_depth,
   const str_vec_t& passive_list) const noexcept;

  /// Setup arrays used throughout `compute`. This includes both arrays that
  /// wrap Cello fields AND temporary arrays used as scratch space.
  ///
  /// @param[in]  block holds data to be processed
  /// @param[out] primitive_map Map of arrays wrapping the Cello Fields holding
  ///     each of the integrable and reconstructable quantity. This also holds
  ///     temporary arrays where the specific form of the passively advected
  ///     scalars will be stored.
  /// @param[out] temp_primitive_map Map for storing the integrable and
  ///     reconstructable quantities at the half timestep. This should have all
  ///     the same entries as primitive_map. However, all arrays in this map
  ///     are temporary.
  /// @param[out] priml_map,primr_map Maps of arrays used to temporarily
  ///     hold the left/right reconstructed face-centered reconstructable and
  ///     integrable quantities. These arrays should have the shape of a
  ///     cell-centered field so that they can be reused for multiple
  ///     dimensions. These have the same keys as primitive_map.
  /// @param[out] pressure_l,pressure_r Arrays used to temporarily store the
  ///     left/right pressure values. The shape of these arrays should be the
  ///     same as in priml_map,primr_map. Note, for adiabatic equations of
  ///     state pressure_l (pressure_r) is a shallow copy of an array in
  ///     priml_map (primr_map).
  /// @param[out] xflux_map, yflux_map, zflux_map Maps of temporary arrays that
  ///     are used to store the x, y, and z fluxes. A given map of arrays will
  ///     hold values along at the face-centers along the direction of the
  ///     fluxes. Note, if a cell-centered field holds `N` elements along
  ///     `dim`, then this should only hold `N-1` elements along `dim`.
  /// @param[out] dUcons_map Map of temporary arrays used to accumulate the
  ///     changes to the conserved forms of the integrable quantities and
  ///     passively advected scalars. If CT is used, this grouping won't have
  ///     space to store changes in the magnetic fields (that update is handled
  ///     separately).
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

  /// Pointer to the BfieldMethod handler
  EnzoBfieldMethod *bfield_method_;

  /// Names of the integrable fields (only includes the field names for
  /// actively advected quantities). These also serve as the keys to the
  /// mappings of arrays used in the calculation
  std::vector<std::string> integrable_field_list_;
  /// Names of the reconstructable primitive fields. These also serve as the
  /// keys to the mappings of arrays used in the calculation
  std::vector<std::string> reconstructable_field_list_;

  /// Lazy initializer of the list of fields holding passive scalars
  EnzoLazyPassiveScalarFieldList lazy_passive_list_;
};

#endif /* ENZO_ENZO_METHOD_VLCT_HPP */
