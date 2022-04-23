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
///    EnzoIntegrationQuanUpdate: Handles updates to integration quantites
///
/// This class can be run with and without magnetic fields. When run with
/// magnetic fields, it makes use of the following component:
///    EnzoBfieldMethod:  Performs constrained transport
///
/// Some notes on implementation
///    - this Method tracks specific total energy (referred to as total_energy)
///    - We categorize quantities as integration quantities quantities and
///      primitives. Nearly every field overlaps between the two cases. Each
///      integration quantity is either a conserved quantity or a conserved
///      quanity divided by density (i.e. specific) Examples of the integration
///      quantities and primitives for adiabatic, ideal gas (without dual
///      energy formalism) include:
///        - density         (both integration and primitive)
///        - velocity        (both integration and primitive)
///        - pressure        (just primitive)
///        - total_energy    (just integration)
///      When supporting magnetic fields, there is also:
///        - bfield          (both integration and primitive)
///      When using the dual energy formalism there is also:
///        - internal_energy (just integration)
///
///    EnzoEFltArrayMap Objects
///    ------------------------
///    The implementation relies upon passing around instances of
///    EnzoEFltArrayMap between the methods of the various components. This
///    class implements a map/dictionary that holds instances of EFlt3DArray.
///    All arrays in a given map are assumed to have the same shape.
///
///    Currently, 9 different maps are used:
///        1. integration_map: Map of arrays wrapping the Cello Fields holding
///           each of the integration quantities. This includes each of the
///           passive scalars (as densities).
///        2. temp_integration_map: Map of arrays used to hold temporary values
///           of for each of the quantities in integration_map. These arrays
///           are used to store the estimated values at the partial timestep.
///        3. primitive_map: Map of arrays used to temporarily store the
///           primitive quantities 
///        4. priml_map: holds left reconstructed primitive quantities (has the
///           same keys as primitive_map)
///        5. primr_map: holds right reconstructed primitive quantities (has the
///           same keys as primitive_map)
///        6. xflux_map: holds fluxes in the x-direction
///        7. yflux_map: holds fluxes in the y-direction
///        8. zflux_map: holds fluxes in the z-direction
///        9. dUcons_map: holds arrays which are used to accumulate the total
///           change in the conserved versions of all integration quantities
///           (includes both flux divergence and source terms)
///    Note:
///        - All arrays in maps 1, 2, 3, and 8 have the same shapes as
///          cell-centered Cello Fields
///        - When they're allocated, all arrays in priml_map and primr_map
///          technically have the shape of a cell-centered field (to allow them
///          to be reused while computing the flux along each dimesnion). While
///          computing the fluxes, the arrays are sliced so that they have the
///          same shape as xflux_map, yflux_map, or zflux_map (depending on the
///          context).
///        - For the purposes of these enumerated maps, we assume that the
///          length of a face-centered array along the dimension with
///          face-centering is 1 less than that of a cell-centered array

#ifndef ENZO_ENZO_METHOD_VLCT_HPP
#define ENZO_ENZO_METHOD_VLCT_HPP

struct EnzoVlctScratchSpace; // defined at the end of this header file

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
		    double dual_energy_formalism_eta,
		    bool store_fluxes_for_corrections);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodMHDVlct);

  /// Charm++ PUP::able migration constructor
  EnzoMethodMHDVlct (CkMigrateMessage *m)
    : Method (m),
      eos_(nullptr),
      half_dt_recon_(nullptr),
      full_dt_recon_(nullptr),
      riemann_solver_(nullptr),
      integration_quan_updater_(nullptr),
      scratch_space_(nullptr),
      mhd_choice_(bfield_choice::no_bfield),
      bfield_method_(nullptr),
      integration_field_list_(),
      primitive_field_list_(),
      lazy_passive_list_(),
      store_fluxes_for_corrections_(false)
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
  virtual double timestep ( Block * block) throw();

protected: // methods

  /// returns the bfield_choice enum that matches the input string
  bfield_choice parse_bfield_choice_(std::string choice) const noexcept;

  /// Checks that the mesh size is sufficiently large to handle the given ghost
  /// depth and confirms that the ghost depth is consistent with the
  /// requirements of the reconstructors
  ///
  /// @param[in] block used to determine the current mesh size and ghost depth
  void check_mesh_and_ghost_size_(Block *block) const noexcept;

  /// Constructs a map containing the field data for each integration quantity
  /// This includes all passively advected scalars (as densities) included in
  /// passive_list
  EnzoEFltArrayMap get_integration_map_(Block * block,
                                        const str_vec_t *passive_list)
    const noexcept;

  /// Computes the fluxes along a given dimension, `dim`, and accumulate the
  /// changes to the integration quantities in `dUcons_map`
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
  /// @param[in]     primitive_map Map of arrays holding cell-centered
  ///     primitive quantities that are to be reconstructed (This includes
  ///     specific passive scalars).
  /// @param[in]     priml_map,primr_map Maps of arrays used to temporarily
  ///     hold the left/right reconstructed face-centered primitives. These
  ///     arrays should have the shape as flux_map
  /// @param[in]     flux_map Holds arrays where the calculated fluxes
  ///     will be stored. The arrays should be face-centered along `dim`.
  ///     If a cell-centered field holds `N` elements along `dim`, then this
  ///     should only hold `N-1` elements along `dim`.
  /// @param[in,out] dUcons_map Map of arrays where the changes to the
  ///     integration quantities are accumulated. If constrained transport is
  ///     being used, this won't include arrays for the magnetic fields.
  /// @param[in]     interface_velocity_arr_ptr Pointer to an array to
  ///     temporarily hold the computed component of the velocity at the cell
  ///     interfaces along `dim`. If a cell-centered field holds `N` elements
  ///     along `dim`, then this is only used to store `N-1` elements. This
  ///     quantity is used to compute the internal energy source term (needed
  ///     under the dual energy formalism). If the value is `nullptr`, then the
  ///     interface velocity is not stored in the array.
  /// @param[in]     reconstructor the instance of EnzoReconstructor to use to
  ///     update reconstruct the face-centered primitives
  /// @param[in,out] bfield_method When using running with bfield handling, this
  ///     is a pointer to an instance of EnzoBfieldMethod. During the function
  ///     call, the internal state is updated. If not handling bfields, this
  ///     should be a `nullptr`.
  /// @param[in]     stale_depth indicates the current stale depth (before
  ///     performing reconstruction)
  /// @param[in]     passive_list A list of keys for passively advected scalars.
  void compute_flux_
  (const int dim, const double cur_dt, const enzo_float cell_width,
   EnzoEFltArrayMap &primitive_map,
   EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
   EnzoEFltArrayMap &flux_map, EnzoEFltArrayMap &dUcons_map,
   const EFlt3DArray* const interface_velocity_arr_ptr,
   EnzoReconstructor &reconstructor, EnzoBfieldMethod *bfield_method,
   const int stale_depth, const str_vec_t& passive_list) const noexcept;

  /// Saves the fluxes (for a given dimension, `dim`), computed at the faces
  /// between the active and the ghost zones to `block->data()->flux_data()`
  /// for later use in flux corrections.
  ///
  /// This function technically saves the value of the flux multiplied by dt
  /// and divided by the cell_width along dimension `dim`.
  ///
  /// @note
  /// This does not currently support flux-correction equivalents for magnetic
  /// fields.
  ///
  /// @param[out] block holds the flux_data object where the fluxes are saved.
  /// @param[in]  flux_map contains the fluxes along that are to be saved.
  /// @param[in]  dim indicates the dimension that the fluxes in `flux_map`
  ///     were computed along. Values of 0, 1, and 2 correspond to the x, y,
  ///     and z directions, respectively.
  /// @param[in]  cell_width is the width of a cell along the dimension `dim`
  /// @param[in]  dt is the value of the current timestep
  void save_fluxes_for_corrections_
  (Block * block, const EnzoEFltArrayMap &flux_map, int dim, double cell_width,
   double dt) const noexcept;

  /// Returns a pointer to the scratch space struct. If the scratch space has
  /// not already been allocated, it will be allocated now.
  ///
  /// This method should be called in EnzoMethodMHDVlct::compute. If it get's
  /// called before the constructors for all methods and initializers are
  /// executed, there's a chance that passive_list may subsequently change.
  ///
  /// @param[in] field_shape Gives the shape, including ghost-zones, of a hydro
  ///     cell-centered field, ordered as (mz,my,mx)
  /// @param[in] passive_list A list of keys for passively advected scalars.
  EnzoVlctScratchSpace* get_scratch_ptr_(const std::array<int,3>& field_shape,
					 const str_vec_t& passive_list)
    noexcept;

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
  /// Pointer to the integration quantity updater
  EnzoIntegrationQuanUpdate *integration_quan_updater_;
  /// Pointer to lazily initialized struct holding scratch-space
  EnzoVlctScratchSpace *scratch_space_;

  /// Indicates how magnetic fields are handled
  bfield_choice mhd_choice_;

  /// Pointer to the BfieldMethod handler
  EnzoBfieldMethod *bfield_method_;

  /// Names of the integration fields (only includes the field names for
  /// actively advected quantities). These also serve as the keys to the
  /// mappings of arrays used in the calculation
  std::vector<std::string> integration_field_list_;
  /// Names of the primitive fields (this should exclude passively advected
  /// scalars). These also serve as the keys to the mappings of arrays used in
  /// the calculation
  std::vector<std::string> primitive_field_list_;

  /// Lazy initializer of the list of fields holding passive scalars
  EnzoLazyPassiveScalarFieldList lazy_passive_list_;

  /// Indicates whether fluxes should be stored for flux corrections
  bool store_fluxes_for_corrections_;
};



struct EnzoVlctScratchSpace{

  /// @class    EnzoVlctScratchSpace
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Holds scratch space arrays for EnzoMethodMHDVlct

public:

  /// Create a new EnzoVlctScratchSpace object
  ///
  /// @param[in] shape Gives the shape, including ghost-zones, of a hydro
  ///     cell-centered field, ordered as (mz,my,mx)
  /// @param[in] integration_key_list List of keys (in the desired order) that
  ///     are associated with each actively-advected cell-centered integration
  ///     quantity. These are used to initialize ``temp_integration_map`` and
  ///     the flux arraymaps
  /// @param[in] primitive_key_list List of keys (in the desired order) that
  ///     are associated with each (non-passively advected) cell-centered
  ///     quantity. These are used to initialize ``primitive_map``,
  ///     ``priml_map`` and ``primr_map``.
  /// @param[in] integ_updater_keys is a list of keys used to initialize
  ///     ``dUcons_map``. This should be returned by the ``integration_keys``
  ///     method of ``EnzoIntegrationQuanUpdate``.
  /// @param[in] passive_list The list of keys for the passively advected
  ///     scalars that should be included in each arraymap.
  /// @param[in] dual_energy Indicates whether the dual energy formalism is in
  ///     use (which specifies if relevant scratch-space should be allocated).
  EnzoVlctScratchSpace(const std::array<int,3>& shape,
                       const str_vec_t& integration_key_list,
                       const str_vec_t& primitive_key_list,
                       const str_vec_t& integ_updater_keys,
                       const str_vec_t& passive_list,
		       bool dual_energy) noexcept
    : interface_vel_arr((dual_energy) ? EFlt3DArray(shape[0],shape[1],shape[2])
			: EFlt3DArray())
  {
    // define function to setup the arraymaps
    auto setup = [&shape, &passive_list](const std::string& name,
                                         const std::array<int,3>& centering,
                                         const str_vec_t& main_keys){
      str_vec_t all_keys(main_keys); // deepcopy of main_keys
      all_keys.insert(all_keys.end(), passive_list.begin(), passive_list.end());
      std::array<int,3> cur_shape(shape); // deepcopy of shape
      for (std::size_t i = 0; i<3; i++){ cur_shape[i] += centering[i]; }
      return EnzoEFltArrayMap(name, all_keys, cur_shape);
    };

    temp_integration_map = setup("temp_integration", {0,0,0},
                                 integration_key_list);
    xflux_map = setup("xflux", { 0, 0,-1}, integration_key_list);
    yflux_map = setup("yflux", { 0,-1, 0}, integration_key_list);
    zflux_map = setup("zflux", {-1, 0, 0}, integration_key_list);
    dUcons_map = setup("dUcons", {0,0,0}, integ_updater_keys);
    primitive_map = setup("primitive", {0,0,0}, primitive_key_list);
    priml_map = setup("priml", {0,0,0}, primitive_key_list);
    primr_map = setup("primr", {0,0,0}, primitive_key_list);
  }

public: // attributes
  /// Array used to store interface velocity values that are computed by the
  /// Riemann Solver(to use in the calculation of the internal energy source
  /// term). If not using the dual energy formalism, this isn't doesn't
  /// allocate memory.
  const CelloArray<enzo_float,3> interface_vel_arr;

  /// Map for storing the integration quantities at the half timestep
  EnzoEFltArrayMap temp_integration_map;

  /// Map of arrays used to temporarily store the cell-centered primitive
  /// quantities that are subsequently reconstructed. This includes arrays for
  /// storing the specific form of each of the passively advected scalars.
  EnzoEFltArrayMap primitive_map;

  /// Maps of arrays used to temporarily hold the left/right reconstructed
  /// face-centered primitive quantities. These have the same keys as
  /// primitive_map. The arrays should have the same shapes as the arrays in
  /// primitive_map (i.e. a shape of cell-centered field), so that they can be
  /// reused for multiple dimensions. As necessary, we take slices so that the
  /// contained arrays are centered along:
  ///   - z and have shape (mz-1,  my,  mx)
  ///   - y and have shape (  mz,my-1,  mx)
  ///   - x and have shape (  mz,  my,mx-1)
  /// where (mz,my,mx) is the shape of an cell-centered array.
  EnzoEFltArrayMap priml_map, primr_map;

  /// Maps of arrays that are used to store the x, y, and z fluxes. If a
  /// cell-centered array has shape (mz,my,mx), then these respectively have
  /// shapes of (mz,my,mx-1), (mz,my-1,mx), and (mz-1,my,mx).
  EnzoEFltArrayMap xflux_map, yflux_map, zflux_map;

  /// Map of temporary arrays used to accumulate the changes to the conserved
  /// forms of the integration quantities and passively advected scalars. If CT
  /// is used, this map won't hold arrays for accumulating changes to the
  /// magnetic fields (that update is handled separately).
  EnzoEFltArrayMap dUcons_map;
};

#endif /* ENZO_ENZO_METHOD_VLCT_HPP */
