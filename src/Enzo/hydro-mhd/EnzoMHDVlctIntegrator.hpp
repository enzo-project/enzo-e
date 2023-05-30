// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMHDVlctIntegrator.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 28 2023
/// @brief    [\ref Enzo] Declaration of EnzoMHDVlctIntegrator

#ifndef ENZO_ENZO_MHD_VLCT_INTEGRATOR_HPP
#define ENZO_ENZO_MHD_VLCT_INTEGRATOR_HPP

class EnzoMHDVlctIntegrator {

  /// @class    EnzoMHDVlctIntegrator
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Solve equations used in for VL + CT MHD method
  ///
  /// This class is intended to encapsulate all algorithmic details about the
  /// the Hydro/MHD solver, other than time-integration details.
  ///
  /// A given Hydro/MHD solver updates hydro-quantites over a single timestep
  /// in one or more stages. The number of stages is determined by the
  /// time-integration choice. For example:
  /// - a simple forward Euler method involves a single stage.
  /// - the van Leer, predictor-corrector method, (which was the original
  ///   choice of solver involves 2 stage).
  /// - Runge-Kutta methods typically involve 2 or more stages.
  ///
  /// The precise details of how such stages are linked together are left to
  /// the disgression of the Method object that drives the integrator. The
  /// Method object is also responsible for taking the Enzo-E fields and
  /// putting them into the correct format for the integrator.
  ///
  /// This class is responsible for tying together all of the algorithmic steps
  /// to update the fluid quantities over an individual stage. It should
  /// basically know nothing about the time integration strategy. There are 2
  /// exceptions to this rule are:
  ///    - it needs to which reconstruction algorithm to use during a given
  ///      stage
  ///    - it may need to know which source terms to skip/compute during a
  ///      given stage.
  /// This class can also provides queryable methods to determine requirements/
  /// assumptions (e.g. one might use them to determine which fields need to be
  /// defined)
  ///
  /// A couple of other design principles for this class include:
  ///    - the internals of this class should be immutable (once an instance is
  ///      constructed, it can't be mutated. In other words, this object
  ///      doesn't carry any mutatable state.
  ///    - this class and its components should know as little as possible
  ///      about the rest of Enzo-E (e.g. it shouldn't directly know anything
  ///      about other Methods, the Block class, the Field class, etc.). The
  ///      main exception is that is allowed to interface with physics objects.
  ///
  /// @note
  /// Currently, the main blemish in the design is that the BfieldMethod class
  /// hierarchy carries a lot of state (it's effectively a state machiner) and
  /// it requires knowledge about the Enzo-E Block class. For that reason,
  /// instances of the BfieldMethod class hierarchy are not currently stored
  /// internally. As a short-term solution, this class provides a method for
  /// constructing the appropriate instance of that class and it must be
  /// explicitly passed an instance of that class for certain operations. A
  /// longer term solution is being developed.
  ///
  /// @note
  /// Going forward, it probably makes sense to move away from implementing the
  /// internal components as a full class hierarchy.

public:

  EnzoMHDVlctIntegrator
  (EnzoRiemann* riemann_solver,
   EnzoReconstructor* half_dt_recon,
   EnzoReconstructor* full_dt_recon,
   EnzoIntegrationQuanUpdate *integration_quan_updater)
    : riemann_solver_(riemann_solver),
      reconstructors_(),
      integration_quan_updater_(integration_quan_updater)
  {
    std::unique_ptr<EnzoReconstructor> half_dt(half_dt_recon);
    std::unique_ptr<EnzoReconstructor> full_dt(full_dt_recon);
    reconstructors_.push_back(std::move(half_dt));
    reconstructors_.push_back(std::move(full_dt));
  }

  /// Delete EnzoMHDVlctIntegrator object
  ~EnzoMHDVlctIntegrator();

  void compute_update_stage
  (EnzoEFltArrayMap tstep_begin_integration_map,
   EnzoEFltArrayMap cur_stage_integration_map,
   EnzoEFltArrayMap out_integration_map,
   EnzoEFltArrayMap primitive_map,
   EnzoEFltArrayMap priml_map,
   EnzoEFltArrayMap primr_map,
   std::array<EnzoEFltArrayMap, 3> flux_maps_xyz,
   EnzoEFltArrayMap dUcons_map,
   const EnzoEFltArrayMap accel_map,
   EFlt3DArray interface_vel_arr,
   const str_vec_t& passive_list,
   EnzoBfieldMethod* bfield_method_,
   unsigned short step_index,
   double cur_dt,
   int stale_depth,
   const std::array<enzo_float,3> cell_widths_xyz) const noexcept;

  /// return the amount that stale_depth increases during a given stage
  ///
  /// @note
  /// Currently, we just consider contributions from reconstruction, but we
  /// could plausibly encounter other sources in the future
  int staling_from_stage(int stage_index) const noexcept {
    check_valid_stage_index_(stage_index);
    return reconstructors_[stage_index]->total_staling_rate();
  }

protected:

  void check_valid_stage_index_(int stage_index) const noexcept {
    int max_size = static_cast<int>(reconstructors_.size());

    ASSERT1("EnzoMHDVlctIntegrator::staling_from_stage",
            "stage_index must satisfy 0 <= stage_index < %d\n",
            max_size, (stage_index >= 0) && (stage_index < max_size));
  }

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
  /// computed in this function AND can be dimensionally split. Other source
  /// terms should be added to the `compute_source_terms_` method.
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

  /// Computes source terms and accumulate the changes to the integration
  /// quantities in `dUcons_map``dU_cons` accordingly.
  ///
  /// At this time, the source terms handled by this method are fairly limited,
  /// but we expect them to grow over time. Note, that the `compute_flux_`
  /// method computes a subset of source terms that can be dimensionally-split
  /// and explicitly depend on relevant dimensional quantities computed by that
  /// method. See the description of that method for more details.
  ///
  /// @param[in]     cur_dt The current timestep.
  /// @param[in]     full_timestep Indicates whether this method is being
  ///     called during the full timestep.
  /// @param[in]     orig_integration_map Map of arrays holding integration
  ///     quantities from the start of the current timestep (this argument
  ///     should be unchanged when calling this method for the partial and then
  ///     full timestep). This nominally includes passive scalars in conserved
  ///     form.
  /// @param[in]     primitive_map Map of arrays holding the current values of
  ///     the primitives (this SHOULD change when calling this method for the
  ///     partial and then full timestep). This nominally includes the
  ///     specific passive scalars
  /// @param[in]     accel_map Map that optionally holds arrays corresponding
  ///     to thr components of the acceleration vector field. This should
  ///     either hold no entries or 3 entries associated with the keys:
  ///     `"acceleration_x"`, `"acceleration_y"`, and `"acceleration_z"`.
  /// @param[in,out] dU_cons Map of arrays where the changes to the
  ///     integration quantities are accumulated.
  /// @param[in]     stale_depth indicates the current stale depth (before
  ///     performing reconstruction)
  ///
  /// @note
  /// The interface of this method will almost certainly need to be updated as
  /// additional source terms get introduced.
  void compute_source_terms_
  (const double cur_dt, const bool full_timestep,
   const EnzoEFltArrayMap &orig_integration_map,
   const EnzoEFltArrayMap &primitive_map,
   const EnzoEFltArrayMap &accel_map,
   EnzoEFltArrayMap &dU_cons, const int stale_depth) const noexcept;

public:
  /// Pointer to the Riemann solver
  EnzoRiemann *riemann_solver_;

  // vector of pointers to reconstructors (the different choices correspond to
  // different stages)
  std::vector<std::unique_ptr<EnzoReconstructor>> reconstructors_;

  /// Pointer to the integration quantity updater
  EnzoIntegrationQuanUpdate *integration_quan_updater_;

};





#endif /* ENZO_ENZO_MHD_VLCT_INTEGRATOR_HPP */
