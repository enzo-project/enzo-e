// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannImpl.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 16 2019
/// @brief    [\ref Enzo] Implementation of the RiemannImpl which is a class
/// template that can be specialized to implement various Riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_IMPL_HPP
#define ENZO_ENZO_RIEMANN_IMPL_HPP

#include <pup_stl.h>
#include <cstdint> // used to check that static methods are defined

//----------------------------------------------------------------------

/// @typedef calc_riemann_fluxes_signature
/// @brief   This is the expected function signature that we expect an
///          `ImplStruct` (a class used to implement a Riemann Solver by
///          specializing `EnzoRiemannImpl`) to have for its static public
///          function, `calc_riemann_fluxes`.
///
/// Within the declaration of an `ImplStruct` we expect the following function:
/// @code
///     static void calc_riemann_fluxes
///       (const enzo_float flux_l[], const enzo_float flux_r[],
///        const enzo_float prim_l[], const enzo_float prim_r[],
///        const enzo_float cons_l[], const enzo_float cons_r[],
///        const enzo_float pressure_l, const enzo_float pressure_r,
///        const EnzoAdvectionFieldLUT lut, const int n_keys,
///        const bool barotropic_eos, const enzo_float gamma,
///        const enzo_float isothermal_cs, const bool dual_energy_formalism,
///        const int iz, const int iy, const int ix,
///        EFlt3DArray flux_arrays[], enzo_float scratch_space[],
///        enzo_float &vi_bar);
/// @endcode
/// This function should computes the Riemann Flux at a given cell interface.
/// `flux_l`, `flux_r`, `prim_l`, `prim_r`, `cons_l`, and `cons_r` store the
/// left and right interface values for fluxes, primitive quantities, and
/// conserved quantities. We note that a given primitive quantities can be
/// either specific or conserved while a conserved quantity is guaranteed to be
/// in conserved form.
///
/// The left and right reconstructed pressure values are passed as `pressure_l`
/// and `pressure_r`.
///
/// The value passed to `lut` maps quantity names to indices in each of the
/// aforementioned arrays and `n_keys` indicates the number of elements held in
/// each array.
///
/// The calculated Riemann flux for a quantity stored at index `j` of the above
/// arrays should be stored at `flux_arrays[j](iz,iy,ix)`.
///
/// `scratch_space` serves as a place to temporarily save quantites during the
/// calculation (the size which should be specified by the `ImplStruct`'s
/// `scratch_space_length` static function).
///
/// `barotropic_eos` indicates whether the fluid equation of state is
/// barotropic. If `true`, then `isothermal_cs` is expected to be non-zero and
/// if `false`, then `gamma` is expected to be positive.
///
/// `dual_energy_formalism` indicates whether to compute fluxes for the
/// internal energy.  It is expected to be `false` when `barotropic_eos==false`
///
/// `vi_bar` is the estimate of `velocity_i` at the cell-interface to be used
/// in the internal energy source term (for the dual energy formalism)
typedef void (*calc_riemann_fluxes_signature)
  (const enzo_float[], const enzo_float[], const enzo_float[],
   const enzo_float[], const enzo_float[], const enzo_float[],
   const enzo_float,   const enzo_float, const EnzoAdvectionFieldLUT,
   const int, const bool, const enzo_float, const enzo_float, const bool,
   const int, const int, const int, EFlt3DArray[], enzo_float[], enzo_float&);

//----------------------------------------------------------------------

/// @typedef scratch_space_length_signature
/// @brief   This is the expected function signature that we expect an
///          `ImplStruct` (a class used to implement a Riemann Solver by
///          specializing `EnzoRiemannImpl`) to have for its static public
///          function, `scratch_space_length`.
///
/// Within the declaration of an `ImplStruct` we expect the following function:
/// @code
///     static int scratch_space_length(const int n_keys);
/// @endcode
/// Given `n_keys`, the number of actively advected quantities for which the
/// Riemann Flux must be computed for, this function returns the length of the
/// `scratch_space` array that is expected to be passed to `ImplStruct`'s
/// `calc_riemann_fluxes` static public method.
typedef int (*scratch_space_length_signature)(const int);


//----------------------------------------------------------------------

/// @def      DEFINE_HAS_SIGNATURE_STATIC
/// @brief    Macro that is used to define a struct that is used to check if a
///           given class definition has a public static function with a
///           specified name and function signature
///
/// @param traitsName The name of the struct that will be defined
/// @param func_name  Name of the static function to check for
/// @param signature  The signature that the function should have (this can be
///     a function pointer)
///
/// This macro comes from https://stackoverflow.com/a/23133787 and is used
/// create structs to help check if the template argument of EnzoRiemannImpl
/// has the correct static methods defined.
#define DEFINE_HAS_SIGNATURE_STATIC(traitsName, func_name, signature)       \
  template <typename U>                                                     \
  class traitsName                                                          \
  {                                                                         \
  private:                                                                  \
    template<typename T, T> struct helper;                                  \
    template<typename T>                                                    \
    static std::uint8_t check(helper<signature, &func_name>*);              \
    template<typename T> static std::uint16_t check(...);                   \
  public:                                                                   \
    static                                                                  \
    constexpr bool value = sizeof(check<U>(0)) == sizeof(std::uint8_t);	    \
  }


// Define the has_calc_riemann_fluxes and has_scratch_space_length to help
// check that classes used to implement riemann solvers by specializing
// EnzoRiemannImpl has the appropriate signatures. We explicitly check this to
// provide more useful debugging messages.
DEFINE_HAS_SIGNATURE_STATIC(has_calc_riemann_fluxes, T::calc_riemann_fluxes,
		            calc_riemann_fluxes_signature);
DEFINE_HAS_SIGNATURE_STATIC(has_scratch_space_length, T::scratch_space_length,
		            scratch_space_length_signature);

template <class ImplStruct>
class EnzoRiemannImpl : public EnzoRiemann
{
  /// @class    EnzoRiemannImpl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides implementation of approximate Riemann
  ///           Solvers
  ///
  /// @tparam ImplStruct The struct used to specialize EnzoRiemannImpl. This
  ///     struct must have static public methods called `calc_riemann_fluxes`
  ///     and `scratch_space_length`. See the documentation for
  ///     `calc_riemann_fluxes_signature` and `scratch_space_length_signature`
  ///     for details about the expected signature and role of each function
  ///
  /// EnzoRiemannImpl factors out the repeated code between different
  /// approximate Riemann Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF &
  /// Roe solvers).
  ///
  /// To define a new RiemannSolver using `EnzoRiemann`:
  ///   1. Define a new `ImplStruct` (e.g. `HLLDImpl`).
  ///   2. It might be useful to define an alias name for the specialization of
  ///      `EnzoRiemannImpl` that uses the new `ImplStruct`
  ///      (e.g. `using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;`).
  ///   3. Then add the particlular specialization to enzo.CI (e.g. add the
  ///      line: `PUPable EnzoRiemannImpl<HLLDImpl>;`)
  ///   4. Update `EnzoRiemann::construct_riemann` to construct the Riemann
  ///      Solver when the correct name is specified.
  ///   5. Update the documentation with the name of the newly available
  ///      RiemannSolver
  ///
  /// @note We could simplify the description of `ImplStruct` by making them
  ///     all functors and having them allocate their scratch space in their
  ///     constructors.


  // The following assertions are first checked when the appropriate .def.h
  // (generated by charm) is included
  static_assert(has_calc_riemann_fluxes<ImplStruct>::value,
		"The static calc_riemann_fluxes method of ImplStruct "
                "doesn't have the right function signature");
  static_assert(has_scratch_space_length<ImplStruct>::value,
		"The static scratch_space_length method of ImplStruct "
                "doesn't have the right function signature");

public: // interface

  /// Constructor
  ///
  /// @param integrable_groups A vector of integrable quantities (listed as
  ///     advected quantities in FIELD_TABLE). These are used as group names in
  ///     the Grouping objects that store field names. In effect this is used
  ///     to register the quantities operated on by the Riemann Solver
  /// @param passive_groups A vector with the names of the groups of passively
  ///     advected scalars that may be included. (If a group is listed here but
  ///     the Grouping object doesn't actually provide any fields in the group,
  ///     no problems are caused
  EnzoRiemannImpl(std::vector<std::string> integrable_groups,
		  std::vector<std::string> passive_groups);

  /// Virtual destructor
  virtual ~EnzoRiemannImpl(){ };

  /// CHARM++ PUP::able declaration
  PUPable_decl_template(EnzoRiemannImpl<ImplStruct>);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemannImpl (CkMigrateMessage *m)
    : EnzoRiemann(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void solve (Block *block, Grouping &priml_group, Grouping &primr_group,
	      std::string pressure_name_l, std::string pressure_name_r,
	      Grouping &flux_group, int dim, EnzoEquationOfState *eos,
	      int stale_depth, std::string interface_velocity_name) const;

protected : //methods
  
  /// Computes the fluxes for the passively advected quantites.
  void solve_passive_advection_(Block* block, Grouping &priml_group,
				Grouping &primr_group, Grouping &flux_group,
				EFlt3DArray &density_flux, int dim,
				int stale_depth) const throw();

  /// Computes the conserved counterpart for every integrable primitive.
  /// Integrable primitives are categorized as conserved, specific, and other
  ///
  /// quantities classified as conserved are copied, quantities classified as
  /// primitive are multiplied by density, and for simplicity, quantities
  /// classified as other are copied (There are no obvious cases where there
  /// should ever be a quanitity classified as other
  void compute_cons_(const enzo_float prim[],
		     enzo_float cons[]) const throw();

  /// computes fluxes for the basic mhd conserved quantities - density,
  /// momentum, energy, magnetic fields  
  void basic_mhd_fluxes_(const enzo_float prim[],
			 const enzo_float cons[],
			 const enzo_float pressure,
			 enzo_float fluxes[]) const throw();

  /// computes extra fluxes for additional, optional physical quantities not
  /// handled in `basic_mhd_fluxes_` (e.g. like internal energy or cosmic rays)
  void extra_fluxes_(const enzo_float prim[], const enzo_float cons[],
		     enzo_float fluxes[],
		     const bool dual_energy_formalism) const throw();

private: //methods

  /// Helper function that simply sets up the lookup table
  void setup_lut_()
  {
    EnzoCenteredFieldRegistry registry;
    lut_ = registry.prepare_advection_lut(integrable_groups_,
					  conserved_start_, conserved_stop_,
					  specific_start_, specific_stop_,
					  other_start_, other_stop_, n_keys_);
    ASSERT("EnzoRiemannImpl::setup_lut_",
	   ("Currently assuming that none of the advected quantites belong to "
	    "the \"other\" category"), other_start_ == other_stop_);
  }

protected: //attributes

  /// Names of the quantities to advect
  std::vector<std::string> integrable_groups_;

  /// struct lookup-table that maps integrable primitive field names to indices
  EnzoAdvectionFieldLUT lut_;

  /// Number of integrable primitive keys (fields) in lut_
  int n_keys_;

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

//----------------------------------------------------------------------

template <class ImplStruct>
EnzoRiemannImpl<ImplStruct>::EnzoRiemannImpl
(std::vector<std::string> integrable_groups,
 std::vector<std::string> passive_groups)
  : EnzoRiemann()
{
  // Quick sanity check - integrable_groups must have density and velocity
  ASSERT("EnzoRiemannImpl","integrable_groups must contain \"density\"",
	 std::find(integrable_groups.begin(), integrable_groups.end(),
		   "density") != integrable_groups.end());
  ASSERT("EnzoRiemannImpl","integrable_groups must contain \"velocity\"",
	 std::find(integrable_groups.begin(), integrable_groups.end(),
		   "velocity") != integrable_groups.end());

  integrable_groups_ = integrable_groups;
  passive_groups_ = passive_groups;

  setup_lut_();
}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::pup (PUP::er &p)
{
  EnzoRiemann::pup(p);

  p|integrable_groups_;
  if (p.isUnpacking()){
    // avoiding PUPing lookup table
    setup_lut_();
  }

  p|passive_groups_;
}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::solve
(Block *block, Grouping &priml_group, Grouping &primr_group,
 std::string pressure_name_l, std::string pressure_name_r,
 Grouping &flux_group, int dim, EnzoEquationOfState *eos, int stale_depth,
 std::string interface_velocity_name) const
{

  const bool barotropic = eos->is_barotropic();
  // if not barotropic then the following doesn't have to be reasonable
  const enzo_float isothermal_cs = eos->get_isothermal_sound_speed();
  // If barotropic, then the following doesn't have to be a reasonable value
  // (this will have to be tweaked when we introduce species that modify gamma)
  const enzo_float gamma = eos->get_gamma();

  // When barotropic equations of state are eventually introduced, all eos
  // dependencies should be moved up here
  ASSERT("EnzoRiemannImpl::solve", "currently no support for barotropic eos",
	 !barotropic);

  const bool dual_energy_formalism = eos->uses_dual_energy_formalism();

  // sanity check:
  ASSERT("EnzoRiemannImpl::solve",
	 ("lut_.internal_energy must have a valid (non-negative) index when "
	  "using the dual energy formalism. Otherwise it must be -1."),
	 ( ( (lut_.internal_energy >= 0) &&  dual_energy_formalism ) ||
	   ( (lut_.internal_energy ==-1) && !dual_energy_formalism ) ));

  // Load arrays:
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  // First, load in the precomputed pressure 
  EFlt3DArray pressure_array_l, pressure_array_r;
  pressure_array_l = array_factory.assigned_center_from_name(pressure_name_l,
							     dim);
  pressure_array_r = array_factory.assigned_center_from_name(pressure_name_r,
							     dim);

  // Second, load field to store interface velocity (if applicable)
  Field field = block->data()->field();
  EFlt3DArray velocity_i_bar_array;
  const bool store_interface_vel =((interface_velocity_name != "") ||
				   (field.is_field(interface_velocity_name)));

  if (store_interface_vel) {
    velocity_i_bar_array =
      array_factory.assigned_center_from_name(interface_velocity_name, dim);
  } else if (dual_energy_formalism) {
    ERROR("EnzoRiemannImpl::solve",
	  "A valid interface_velocity_name must be provided when the dual "
	  "energy formalism is in use");
  }

  EnzoCenteredFieldRegistry registry;
  EFlt3DArray *wl_arrays, *wr_arrays, *flux_arrays;
  wl_arrays = registry.load_array_of_fields(block, lut_, n_keys_,
					    priml_group, dim, stale_depth);
  wr_arrays = registry.load_array_of_fields(block, lut_, n_keys_,
					    primr_group, dim, stale_depth);
  flux_arrays = registry.load_array_of_fields(block, lut_, n_keys_,
					      flux_group, dim, stale_depth);

  // allocate arrays to temporarily hold values at each cell interface
  enzo_float *wl, *wr, *Ul, *Ur, *Fl, *Fr;
  wl = new enzo_float[n_keys_];    wr = new enzo_float[n_keys_];
  Ul = new enzo_float[n_keys_];    Ur = new enzo_float[n_keys_];
  Fl = new enzo_float[n_keys_];    Fr = new enzo_float[n_keys_];

  // prepare optional scratch space to be used by ImplStruct
  enzo_float *scratch_space = NULL;
  // More efficient code might be generated if we were just to pass the maximum
  // possible number of keys as the argument here. Then the size of scratch
  // space could be determined at compile time (allowing scratch_space to be
  // allocated off the stack)
  int scratch_space_length = ImplStruct::scratch_space_length(n_keys_);
  if (scratch_space_length>0){
    scratch_space = new enzo_float[scratch_space_length];
  }


  // For Nearest-Neighbor, we care about interfaces starting at i+1/2. Thanks
  // to use of stale_depth, we can treat every case like this
  for (int iz=0; iz<flux_arrays[0].shape(0); iz++) {
    for (int iy=0; iy<flux_arrays[0].shape(1); iy++) {
      for (int ix=0; ix<flux_arrays[0].shape(2); ix++) {

	// get the fluid fields
	for (int field_ind=0; field_ind<n_keys_; field_ind++){
	  wl[field_ind] = wl_arrays[field_ind](iz,iy,ix);
	  wr[field_ind] = wr_arrays[field_ind](iz,iy,ix);
	}

	// get the left/right pressure
	enzo_float pressure_l = pressure_array_l(iz,iy,ix);
	enzo_float pressure_r = pressure_array_r(iz,iy,ix);

	// get the conserved quantities
	compute_cons_(wl, Ul);
	compute_cons_(wr, Ur);

	// compute the interface fluxes
	basic_mhd_fluxes_(wl, Ul, pressure_l, Fl);
	basic_mhd_fluxes_(wr, Ur, pressure_r, Fr);

	// compute the other (optional) fluxes
	extra_fluxes_(wl, Ul, Fl, dual_energy_formalism);
	extra_fluxes_(wr, Ur, Fr, dual_energy_formalism);

	enzo_float interface_velocity_i;
	// Now compute the Riemann Fluxes
	ImplStruct::calc_riemann_fluxes(Fl, Fr, wl, wr, Ul, Ur,
				        pressure_l, pressure_r, lut_, n_keys_,
					barotropic, gamma, isothermal_cs,
					dual_energy_formalism, iz, iy, ix,
					flux_arrays, scratch_space,
					interface_velocity_i);

	if (store_interface_vel){
	  velocity_i_bar_array(iz,iy,ix) = interface_velocity_i;
	}

	/*
	// If Dedner Fluxes are required, they get handled here
	//   - It would probably be better to handle this separately from the
	//     Riemann Solver since we already precompute L & R conserved
	//     AND it doesn't require wavespeed information.
	if (Dedner_Ch_ != 0){
	  // Not sure this is perfectly transposed
	  flux_arrays["bfield_i"](iz,iy,ix) =
	    (Ul["phi"] + 0.5*(Ur["phi"] - Ul["phi"])
	     - 0.5 * Dedner_Ch_ * (Ur["bfield_i"] - Ul["bfield_i"]));
	  flux_arrays["phi"](iz,iy,ix) =
	    (Ul["bfield_i"] + 0.5*(Ur["bfield_i"] - Ul["bfield_i"])
	     - 0.5 / Dedner_Ch_ * (Ur["phi"] - Ul["phi"]));
	  flux_arrays["phi"](iz,iy,ix) *= (Dedner_Ch_*Dedner_Ch_);
	}
	*/
      }
    }
  }

  solve_passive_advection_(block, priml_group, primr_group, flux_group,
			   flux_arrays[lut_.density], dim, stale_depth);

  delete[] wl; delete[] wr;
  delete[] Ul; delete[] Ur;
  delete[] Fl; delete[] Fr;
  delete[] wl_arrays; delete[] wr_arrays;
  delete[] flux_arrays;
  if (scratch_space != NULL){
    delete[] scratch_space;
  }
}

inline void compute_unity_sum_passive_fluxes_(const enzo_float dens_flux,
					      EFlt3DArray *flux_arrays,
					      EFlt3DArray *reconstructed,
					      const int num_fields,
					      const int iz, const int iy,
					      const int ix) throw()
{
  enzo_float sum = 0.;
  for (int field_ind=0; field_ind<num_fields; field_ind++){
    sum += reconstructed[field_ind](iz,iy,ix);
  }
  for (int field_ind=0; field_ind<num_fields; field_ind++){
    flux_arrays[field_ind](iz,iy,ix) = (reconstructed[field_ind](iz,iy,ix)
					* dens_flux / sum);
  }
}


template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::solve_passive_advection_
(Block* block, Grouping &priml_group, Grouping &primr_group,
 Grouping &flux_group, EFlt3DArray &density_flux, int dim, int stale_depth)
  const throw()
{
  // This was basically transcribed from Enzo
  std::vector<std::string> group_names = this->passive_groups_;

  EnzoFieldArrayFactory array_factory(block, stale_depth);
  EnzoPermutedCoordinates coord(dim);

  // unecessary values are computed for the inside faces of outermost ghost zone
  for (unsigned int group_ind=0; group_ind<group_names.size(); group_ind++){

    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = priml_group.size(group_name);

    if (num_fields == 0){
      continue;
    }

    // If we have a group_name that we have to make sure sums to 1, we should
    // do that now! (The easiest way to facillitate that would be to separate
    // spieces from colours)
    bool unity_sum = false;

    // load array of fields
    EFlt3DArray *wl_arrays = new EFlt3DArray[num_fields];
    EFlt3DArray *wr_arrays = new EFlt3DArray[num_fields];
    EFlt3DArray *flux_arrays = new EFlt3DArray[num_fields];

    for (int field_ind=0; field_ind<num_fields; field_ind++){
      wl_arrays[field_ind] = array_factory.reconstructed_field(priml_group,
							       group_name,
							       field_ind, dim);
      wr_arrays[field_ind] = array_factory.reconstructed_field(primr_group,
							       group_name,
							       field_ind, dim);
      flux_arrays[field_ind] = array_factory.from_grouping(flux_group,
							   group_name,
							   field_ind);
    }

    for (int iz=0; iz<density_flux.shape(0); iz++) {
      for (int iy=0; iy<density_flux.shape(1); iy++) {
	for (int ix=0; ix<density_flux.shape(2); ix++) {

	  enzo_float dens_flux = density_flux(iz,iy,ix);
	  EFlt3DArray *reconstr = (dens_flux>0) ? wl_arrays : wr_arrays;

	  if (unity_sum) {
	    compute_unity_sum_passive_fluxes_(dens_flux, flux_arrays, reconstr,
					      num_fields, iz, iy, ix);
	  } else {
	    for (int field_ind=0; field_ind<num_fields; field_ind++){
	      flux_arrays[field_ind](iz,iy,ix) = (reconstr[field_ind](iz,iy,ix)
						  * dens_flux);
	    }
	  }

	}
      }
    }
    delete[] wl_arrays; delete[] wr_arrays; delete[] flux_arrays;
  }
}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::compute_cons_
(const enzo_float prim[], enzo_float cons[]) const throw()
{
  enzo_float density = prim[lut_.density];

  // includes density, bfields, ...
  for (int i= conserved_start_; i<conserved_stop_; i++){
    cons[i] = prim[i];
  }

  // includes velocities, specific total energy, ...
  for (int i= specific_start_; i<specific_stop_; i++){
    cons[i] = density * prim[i];
  }

  ASSERT("EnzoRiemannImpl::compute_cons_",
	 "As of now other_start_ should always be equal to other_stop_",
	 other_start_ == other_stop_);
  //for (int i= other_start_; i<other_stop_; i++){
  //  cons[i] = prim[i];
  //}
}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::basic_mhd_fluxes_
(const enzo_float prim[], const enzo_float cons[], const enzo_float pressure,
 enzo_float fluxes[]) const throw()
{
  // Assumes that MHD is included. Not currently equipped for barotropic eos. 
  // This may be better handled by the EquationOfState

  // when dual-energy formalism is in use, the internal energy fluxes are
  // handeled by EnzoRiemannImpl<ImplStruct>::extra_fluxes_

  enzo_float vi, vj, vk, p, Bi, Bj, Bk, etot, mag_pressure;
  vi = prim[lut_.velocity_i];
  vj = prim[lut_.velocity_j];
  vk = prim[lut_.velocity_k];
  Bi = prim[lut_.bfield_i];
  Bj = prim[lut_.bfield_j];
  Bk = prim[lut_.bfield_k];
  etot = cons[lut_.total_energy];

  p  = pressure;

  mag_pressure = mag_pressure_(prim, lut_);

  // Compute Fluxes
  enzo_float pi = cons[lut_.velocity_i];
  fluxes[lut_.density] = pi;

  // Fluxes for Mx, My, Mz
  fluxes[lut_.velocity_i] = pi*vi - Bi*Bi + p + mag_pressure;
  fluxes[lut_.velocity_j] = pi*vj - Bj*Bi;
  fluxes[lut_.velocity_k] = pi*vk - Bk*Bi;

  // Flux for etot
  fluxes[lut_.total_energy] = ((etot + p + mag_pressure)*vi
			       - (Bi*vi + Bj*vj + Bk*vk)*Bi);

  // Fluxes for Bi,Bj,Bk
  fluxes[lut_.bfield_i] = 0;
  fluxes[lut_.bfield_j] = Bj*vi - Bi*vj;
  fluxes[lut_.bfield_k] = Bk*vi - Bi*vk;
}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::extra_fluxes_
(const enzo_float prim[], const enzo_float cons[], enzo_float fluxes[],
 const bool dual_energy_formalism) const throw()
{
  // as more physics get added (e.g. this will need to be extended)
  if (dual_energy_formalism){
    fluxes[lut_.internal_energy] =
      cons[lut_.internal_energy] * prim[lut_.velocity_i];
  }
}

#endif /* ENZO_ENZO_RIEMANN_IMPL_HPP */
