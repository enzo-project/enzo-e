// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannImpl.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 16 2019
/// @brief    [\ref Enzo] Implementation of the RiemannImpl which is a class
/// template that can be specialized to implement various riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_IMPL_HPP
#define ENZO_ENZO_RIEMANN_IMPL_HPP

#include <pup_stl.h>
#include <cstdint> // used to check that static methods are defined

// EnzoRiemannImpl class factors out the repeated code between different
// approximate Riemann Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF & Roe
// solvers). To facilitate this, EnzoRiemannImpl is a class template.
//
// To implement a specific solver one needs to define an ImplStruct, a struct
// that must have have the following 2 static methods:
//
//  static void calc_riemann_fluxes
//   (const enzo_float flux_l[], const enzo_float flux_r[],
//    const enzo_float prim_l[], const enzo_float prim_r[],
//    const enzo_float cons_l[], const enzo_float cons_r[],
//    const enzo_float pressure_l, const enzo_float pressure_r,
//    const EnzoAdvectionFieldLUT lut, const int n_keys,
//    const bool barotropic_eos, const enzo_float gamma,
//    const enzo_float isothermal_cs,
//    const int iz, const int iy, const int ix,
//    EFlt3DArray flux_arrays[], enzo_float scratch_space[]);
//
//   This function computes the Riemann Flux at a given cell interface. flux_l,
//   flux_r, prim_l, prim_r, cons_l, cons_r, indicate the values the left and
//   right interface values for fluxes, conserved quantities, and primitives.
//   prim_lut and cons_lut map quantity names to indices in these arrays. n_keys
//   indicates the number of conserved quantities. The resulting Riemann flux
//   for conserved quantity at index j get's stored at flux_arrays[j](iz,iy,ix)
//   Finally, scratch_space serves as a place to temporarily save quantites
//   during the calculation. If barotropic_eos is True, then
//   isothermal_cs is expected to be non-zero, while if it's False, then gamma
//   is expected to be positive.
//
//  static int scratch_space_length(const int n_cons_keys);
//
//  This function determines the size of scratch space required given the
//  number of conserved quantities.
//
// To define a new RiemannSolver:
// 1. Define a new ImplStruct (e.g. HLLDImpl).
// 2. It might be useful to define an alias name for the specialization of
//    EnzoRiemannImpl that uses the new ImplStruct
//    (e.g. "using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;").
// 3. Then add the particlular specialization to enzo.CI (e.g. add the line
//    "PUPable EnzoRiemannImpl<HLLDImpl>;")
// 4. Finally, update EnzoRiemann::construct_riemann to construct the Riemann
//    Solver under the appropriate conditions


class EnzoFluxFunctor : public PUP::able
{
  /// @class    EnzoFluxFunctor
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Allows for easy implementation of additional
  ///           non-passively advected (like cosmic rays). Subclasses of this
  ///           object should be passed to RiemannImpl's constructor
  ///
  /// Provides an easy way to add additional fields to RiemannImpl and
  /// derivatives. The resulting functors simply get called to compute Riemann
  /// fluxes
  ///
  /// It's almost certain that these will be too slow due to the indirection.
  /// There are 2 main alternatives:
  ///     - could hardcode the functions into the RiemannImpl or passing some
  ///       variadic template arguments to RiemannSolver::solve (hardcoding
  ///       may be the way to go)
  
public:
  EnzoFluxFunctor() throw()
  {}

  virtual ~EnzoFluxFunctor()
  {}

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoFluxFunctor);

   /// CHARM++ migration constructor for PUP::able
  EnzoFluxFunctor (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  virtual void operator()(const enzo_float prim[], const enzo_float cons[],
			  enzo_float fluxes[], const EnzoRiemannImpl lut)=0;
};

// This should all get moved to separate documentation page on website
//
// Below is a brief summary of EnzoRiemannImpl and a brief description of how
// EnzoAdvectionFieldLUT struct, derived from the FIELD_TABLE XMacro, is used
//   - The Riemann Solver makes use of two sets of arrays.
//
//       1. The first set is made up of arrays of instances of EFlt3DArray.
//          Within a given array, each instance of EFlt3DArray encapsulates the
//          data for a different field. The solver uses:
//            A. arrays of left/right reconstructed integrable primitive fields
//            B. arrays of fields where the calculated Riemann Flux is stored
//
//       2. The second set is made up of arrays of enzo_float. These arrays
//          serve as temporary storage buffers which store quantities for a
//          single cell interface. The solver uses
//            A. arrays of left/right reconstructed integrable primitives
//            B. arrars of left/right fluxes
//            C. arrays of left/right reconstructed conserved quantities
//               (these are computed from the integrable primitives, which may
//                conserved, specific, or other)
//
//     The left and right pressures are also used. These are passed in
//     separately and computed ahead of time with the equations of state
//
//     Each of the above arrays only include quantities that are required for
//     the calculations of wave speeds or have non-trivial flux calculations.
//     Passively advected scalars are not included in these arrays (their flux
//     is computed separately).
//
//   - The general flow of the EnzoRiemannImpl::solve is:
//
//       a For a given cell-interface on a grid, the reconstructed integrator
//         primitives are copied from arrays 1A into the temporary array 2A.
//
//       b The conserved quantities in array 2C are calculated/copied from the
//         quantities in array 2A. (Integrable primitives are classified as
//         conserved, specific, and other. For now, we just copy the quantities
//         in other over to the array of quantities in conserved)
//
//       c The standard MHD fluxes are computed at that location and saved into
//         the arrays left/right fluxes.
//
//       d Optional functions are also applied to compute additional left/right
//         fluxes. Pointers to these functions are specified upon construction
//         of the solver and are used for non-standard fluxes (e.g. cosmic ray
//         energy/fluxes)
//
//       e This step is implemented in static functions of the ImplStruct
//         template argument. The wavespeeds at the current interface is
//         are computed, and the total riemann fluxes are computed at the
//         interface. The calculation of the total fluxes are achieved by
//         iterating over the conserved and flux quantities
//
//   - Use of the following EnzoAdvectionFieldLUT struct:
//
//       - The calculation of fluxes and wave speeds requires random access of
//         specific fields. We also need to be able to iterate over the entries
//         of multiple arrays simultaneously (e.g. to accomplish part e, above)
//         Unlike Enzo, we wanted to avoid statically declaring which indices
//         correspond to which fields (adding additional sets of fields, like
//         internal energy and cosmic ray energy/fluxes becomes harder)
//
//       - We settled on using the field_lut struct as a lookup table. The
//         struct has members named for every quantity listed in FIELD_TABLE
//           - For a SCALAR, the member name directly matches the name in
//             column 1
//           - For a VECTOR, there are 3 members: {name}_i, {name}_j, {name}_k
//             ({name} corresponds to the name appearing in column 1)
//         Each struct contains members named advection-related (related to
//         reconstruction or riemann fluxes) quantities in the table.
//
//       - Example: If we have an array of reconstructed integrable primitives,
//         wl, and an initialized instance of EnzoAdvectionFieldLUT, lut,
//         then wl[prim_lut.density] and wl[prim_lut.total_energy]
//         indicates the entries reserved for the density and (specific)
//         total energy
//
//       - EnzoCenteredFieldRegistry::prepare_advection_lut is used to
//         prepare the lookup table for pre-specified quantities and determines
//         the length of an array large enough to hold fields representing all
//         of the quantities. The function then initializes all members of the
//         struct to be equal to quantities from 0 through 1 less than the
//         length of an arrays and organizes them based on whether the quanitiy
//         is conserved, specific, or other. The function also yields the
//         indices required to just iterate over each category of field. All
//         members of the struct that don't correspond to a specified quantity
//         are set to -1
//     
//       - load_array_of_fields constructs an array of instances of EFlt3DArray
//         that correspond to reconstructed integrable primitive fields OR flux
//         fields. The ordering of the fields is determined by the supplied
//         instance of EnzoAdvectionFieldLUT


//----------------------------------------------------------------------
/// @def      DEFINE_HAS_SIGNATURE
/// @brief    Macro that is used to define a struct (whose name is passed as
///           the traitsName) that is used to check if a given class definition
///           has a static function with a specified name (passed as the
///           func_name argument of this macro) and function signature (passed
///           as the signature arguement)
///
/// This macro comes from https://stackoverflow.com/a/23133787 and is used
/// create structs to help check if the template argument of EnzoRiemannImpl
/// has the correct functions defined.
#define DEFINE_HAS_SIGNATURE(traitsName, func_name, signature)              \
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

DEFINE_HAS_SIGNATURE(has_calc_riemann_fluxes, T::calc_riemann_fluxes,
		     void (*)(const enzo_float[], const enzo_float[],
			      const enzo_float[], const enzo_float[],
			      const enzo_float[], const enzo_float[],
			      const enzo_float, const enzo_float,
			      const field_lut, const int,
			      const bool, const enzo_float, const enzo_float,
			      const int, const int, const int,
			      EFlt3DArray[], enzo_float[]));
DEFINE_HAS_SIGNATURE(has_scratch_space_length, T::scratch_space_length,
		     int (*)(const int));

template <class ImplStruct>
class EnzoRiemannImpl : public EnzoRiemann
{
  /// @class    EnzoRiemannImpl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides implementation of approximate Riemann
  ///           Solvers
  ///
  /// @tparam ImplStruct The struct used to specialize EnzoRiemannImpl


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
  /// @param flux_funcs pointer to an array of instances of EnzoFluxFunctor
  /// @param n_funcs the size of the array of flux functors
  EnzoRiemannImpl(std::vector<std::string> integrable_groups,
		  std::vector<std::string> passive_groups,
		  EnzoFluxFunctor** flux_funcs, int n_funcs);

  /// Virtual destructor
  virtual ~EnzoRiemannImpl();

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
	      int stale_depth);

protected : //methods
  
  /// Computes the fluxes for the passively advected quantites.
  void solve_passive_advection_(Block* block, Grouping &priml_group,
				Grouping &primr_group, Grouping &flux_group,
				EFlt3DArray &density_flux, int dim);

  /// Computes the conserved counterpart for every integrable primitive.
  /// Integrable primitives are categorized as conserved, specific, and other
  ///
  /// quantities classified as conserved are copied, quantities classified as
  /// primitive are multiplied by density, and for simplicity, quantities
  /// classified as other are copied (There are no obvious cases where there
  /// should ever be a quanitity classified as other
  void compute_cons_(const enzo_float prim[], enzo_float cons[]);

  
  /// computes fluxes for the basic mhd conserved quantities - density,
  /// momentum, energy, magnetic fields  
  void basic_mhd_fluxes_(const enzo_float prim[], const enzo_float cons[],
			 const enzo_float pressure, enzo_float fluxes[]);
private: //methods

  /// Helper function that simply sets up the lookup table
  void setup_lut_()
  {
    EnzoCenteredFieldRegistry registry;
    lut_ = registry.prepare_advection_lut(integrable_groups_,
					  conserved_start_, conserved_stop_,
					  specific_start_, specific_stop_,
					  other_start_, other_stop_, n_keys);
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

  /// number of flux functors
  int n_funcs_;

  /// array of pointers to functors used to compute fluxes
  EnzoFluxFunctor** flux_funcs_;

};

//----------------------------------------------------------------------

template <class ImplStruct>
EnzoRiemannImpl<ImplStruct>::EnzoRiemannImpl
(std::vector<std::string> integrable_groups,
 std::vector<std::string> passive_groups,
 EnzoFluxFunctor** flux_funcs, int n_funcs)
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

  n_funcs_ = n_funcs;
  flux_funcs_ = flux_funcs;

  setup_lut_();
}

//----------------------------------------------------------------------

template <class ImplStruct>
EnzoRiemannImpl<ImplStruct>::~EnzoRiemannImpl()
{
  for (int i = 0; i < n_funcs_; i++){
    delete flux_funcs_[i];
  }
  delete[] flux_funcs_;
}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::pup (PUP::er &p)
{
  EnzoRiemann::pup(p);

  p|conditions_;
  if (p.isUnpacking()){
    // avoiding PUPing lookup table
    setup_lut_();
  }

  p|passive_groups_;

  p|n_funcs_;

  if (p.isUnpacking()){
    flux_funcs_ = new EnzoFluxFunctor*[n_funcs_];
  }
  for (int i = 0; i<n_funcs_;i++){
    p|flux_funcs_[i];
  }
}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::solve
(Block *block, Grouping &priml_group, Grouping &primr_group,
 std::string pressure_name_l, std::string pressure_name_r,
 Grouping &flux_group, int dim, EnzoEquationOfState *eos, int stale_depth)
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

  ASSERT("EnzoRiemannImpl::solve",
	 "currently no support for dual energy formalism",
	 !eos->uses_dual_energy_formalism());

  // Let's load in the precomputed pressure 
  EFlt3DArray pressure_array_l, pressure_array_r;
  EnzoFieldArrayFactory array_factory(block, stale_depth);
  pressure_array_l = array_factory.reconstructed_from_name(pressure_name_l,
							   dim);
  pressure_array_r = array_factory.reconstructed_from_name(pressure_name_r,
							   dim);

  EnzoCenteredFieldRegistry registry;
  EFlt3DArray *wl_arrays, *wr_arrays, *flux_arrays;

  wl_arrays = registry.load_array_of_fields(block, lut_, n_keys_,
					    priml_group, dim, stale_depth);
  wr_arrays = registry.load_array_of_fields(block, lut_, n_keys_,
					    primr_group, dim, stale_depth);
  flux_arrays = registry.load_array_of_fields(block, lut_, n_keys_,
					      flux_group, dim, stale_depth);

  enzo_float *wl, *wr, *Ul, *Ur, *Fl, *Fr;
  wl = new enzo_float[n_keys_];    wr = new enzo_float[n_keys_];
  Ul = new enzo_float[n_keys_];    Ur = new enzo_float[n_keys_];
  Fl = new enzo_float[n_keys_];    Fr = new enzo_float[n_keys_];

  // prepare optional scratch space to be used by ImplStruct
  enzo_float *scratch_space = NULL;
  int scratch_space_length = ImplStruct::scratch_space_length(n_keys_);
  if (scratch_space_length>0){
    scratch_space = new enzo_float[scratch_space_length];
  }
  

  // For Nearest-Neighbor, we care about interfaces starting at i+1/2
  // For PLM we only care about interfaces starting at i+3/2.
  //    (left interface value at i+1/2 is set to 0)
  // For consistency, always start at i+1/2

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

	// iterate over the functors
	for (int i = 0; i<n_funcs_; i++){
	  (*(flux_funcs_[i]))(wl, Ul, Fl, prim_lut_, cons_lut_);
	  (*(flux_funcs_[i]))(wr, Ur, Fr, prim_lut_, cons_lut_);
	}


	// Now compute the Riemann Fluxes
	ImplStruct::calc_riemann_fluxes(Fl, Fr, wl, wr, Ul, Ur,
				        pressure_l, pressure_r, lut_, n_keys_,
					barotropic, gamma, isothermal_cs, iz,
					iy, ix, flux_arrays, scratch_space);

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
			   flux_arrays[cons_lut_.density], dim);

  delete[] wl; delete[] wr;
  delete[] Ul; delete[] Ur;
  delete[] Fl; delete[] Fr;
  delete[] wl_arrays; delete[] wr_arrays;
  delete[] ul_arrays; delete[] ur_arrays;
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
					      const int ix)
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
 Grouping &flux_group, EFlt3DArray &density_flux, int dim)
{
  // This was basically transcribed from Enzo-E
  std::vector<std::string> group_names = this->passive_groups_;

  EnzoFieldArrayFactory array_factory(block);
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
      wl_arrays[field_ind] = array_factory.from_grouping(priml_group,group_name,
							 field_ind);
      wr_arrays[field_ind] = array_factory.from_grouping(primr_group,group_name,
							 field_ind);
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
void EnzoRiemannImpl<ImplStruct>::compute_cons_(const enzo_float prim[],
						enzo_float cons[])
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

  // I don't think this should include anything
  for (int i= other_start_; i<other_stop_; i++){
    other[i] = prim[i];
  }

}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::basic_mhd_fluxes_
(const enzo_float prim[], const enzo_float cons[], const enzo_float pressure,
 enzo_float fluxes[])
{
  // This assumes that MHD is included
  // not currently equipped for barotropic OR dual energy formalism
  
  // This may be better handled by the EquationOfState
  enzo_float vi, vj, vk, p, Bi, Bj, Bk, etot, mag_pressure;
  vi = prim[lut_.velocity_i];
  vj = prim[lut_.velocity_j];
  vk = prim[lut_.velocity_k];
  Bi = prim[lut_.bfield_i];
  Bj = prim[lut_.bfield_j];
  Bk = prim[lut_.bfield_k];
  etot = cons[lut_.total_energy];

  p  = pressure;

  mag_pressure = mag_pressure_(prim, prim_lut);

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

#endif /* ENZO_ENZO_RIEMANN_IMPL_HPP */
