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
///        - All arrays in priml_map and primr_map technically have the shape
///          of a cell-centered field, but they are treated as though they have
///          the shape of a face-centered field. This is done so that they can
///          be reused for different axes.
///        - For the purposes of these enumerated maps, we assume that the
///          length of a face-centered array along the dimension with
///          face-centering is 1 less than that of a cell-centered array

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
      integration_quan_updater_(nullptr),
      mhd_choice_(bfield_choice::no_bfield),
      bfield_method_(nullptr),
      integration_field_list_(),
      primitive_field_list_(),
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
  ///     arrays should have the shape of a cell-centered field, but are
  ///     treated as though they have the shape of a that is face-centered
  ///     along `dim` (If a cell-centered field holds `N` elements along `dim`,
  ///     then such a face-centered field should only have `N-1` elements along
  ///     `dim`).
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
   EFlt3DArray *interface_velocity_arr_ptr, EnzoReconstructor &reconstructor,
   EnzoBfieldMethod *bfield_method, const int stale_depth,
   const str_vec_t& passive_list) const noexcept;

  /// Setup arrays used throughout `compute`. This includes both arrays that
  /// wrap Cello fields AND temporary arrays used as scratch space.
  ///
  /// @param[in]  block holds data to be processed
  /// @param[out] integration_map Map of arrays wrapping the Cello Fields
  ///     holding each of the integration quantities. This includes each of the
  ///     passive scalars (as densities).
  /// @param[out] temp_integration_map Map for storing the integration
  ///     quantities at the half timestep. This should have all
  ///     the same entries as primitive_map. However, all arrays in this map
  ///     are temporary.
  /// @param[out] primitive_map Map of arrays used to temporarily store the
  ///     cell-centered primitive quantities that are subsequently
  ///     reconstructed. This includes arrays for storing the specific form of
  ///     each of the passively advected scalars.
  /// @param[out] priml_map,primr_map Maps of arrays used to temporarily
  ///     hold the left/right reconstructed face-centered primitive quantities.
  ///     These arrays should have the shape of a cell-centered field so that
  ///     they can be reused for multiple dimensions. These have the same keys
  ///     as primitive_map.
  /// @param[out] xflux_map, yflux_map, zflux_map Maps of temporary arrays that
  ///     are used to store the x, y, and z fluxes. A given map of arrays will
  ///     hold values along at the face-centers along the direction of the
  ///     fluxes. Note, if a cell-centered field holds `N` elements along
  ///     `dim`, then this should only hold `N-1` elements along `dim`.
  /// @param[out] dUcons_map Map of temporary arrays used to accumulate the
  ///     changes to the conserved forms of the integration quantities and
  ///     passively advected scalars. If CT is used, this grouping won't have
  ///     space to store changes in the magnetic fields (that update is handled
  ///     separately).
  void setup_arrays_
  (Block *block, EnzoEFltArrayMap &integration_map,
   EnzoEFltArrayMap &temp_integration_map, EnzoEFltArrayMap &primitive_map,
   EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
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
  /// Pointer to the integration quantity updater
  EnzoIntegrationQuanUpdate *integration_quan_updater_;

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
};

#endif /* ENZO_ENZO_METHOD_VLCT_HPP */
