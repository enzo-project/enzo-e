// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannImpl.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 16 2019
/// @brief    [\ref Enzo] Implementation of the Riemann Solver abstract base
/// class. This class should be subclassed to implement various riemann solvers.

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
//    const field_lut prim_lut, const field_lut cons_lut, const int n_keys,
//    const bool barotropic_eos, const enzo_float gamma,
//    const enzo_float isothermal_cs, const int iz, const int iy, const int ix,
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


// To allow for easily adding additional fields with non-trivial flux
// calculations (e.g. cosmic rays or possibly internal energy), the
// constructor will accept an array of functors, subclassed from FluxFunctor
// - If the functors are too slow, we could hardcode the functions into the
//   Riemann Solver or we could specify the functors as variadic template
//   arguments (allowing operator() to be inlined within the for loop)
// - this should just be replaced with hard coded function within
//   RiemannSolver::solve or some use of variadic templates (maybe a tuple).
//   Hard-coding may be the way to go 
//
// To add an additional field, add the field to FIELD_TABLE and modify the
// factory method

class FluxFunctor : public PUP::able
{
public:
  FluxFunctor() throw()
  {}

  virtual ~FluxFunctor()
  {}

  /// CHARM++ PUP::able declaration
  PUPable_abstract(FluxFunctor);

   /// CHARM++ migration constructor for PUP::able
  FluxFunctor (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  virtual void operator()(const enzo_float prim[], const enzo_float cons[],
			  enzo_float fluxes[], const field_lut prim_lut,
			  const field_lut cons_lut)=0;
};


// Below is a brief summary of the derivatives of the FIELD_TABLE XMacro that
// are employed by EnzoRiemannImpl
//   - The Riemann Solver makes use of two sets of arrays.
//
//       1. The first set is made up of arrays of instances of EFlt3DArray.
//          Within a given array, each instance of EFlt3DArray encapsulates the
//          data for a different field. The solver uses:
//            A. arrays of left/right reconstructed primitive fields
//            B. arrays of left/right reconstructed conserved fields
//            C. array of fields where the calculated Riemann Flux is stored
//
//       2. The second set is made up of arrays of instances of enzo_float.
//          These arrays serve as temporary storage buffers which store 
//          quantities for a single cell interface. The solver uses
//            A. arrays of left/right reconstructed primitives quantities
//            B. arrays of left/right reconstructed conserved quantities
//            C. arrars of left/right fluxes
//
//     Each of the above arrays only include quantities that are required for
//     the calculations of wave speeds or have non-trivial flux calculations.
//     Passively advected scalars are not included in these arrays (their flux
//     is computed separately).
//
//   - The general code-flow of the RiemannSolver is:
//
//       a For a given cell-interface on a grid, the reconstructed primitive
//         conserved quantites from arrays 1A and 1B into the temporary arrays
//         2A and 2B.
//
//       b The standard MHD fluxes are computed at that location and saved into
//         the arrays left/right fluxes.
//
//       c Optional functions are also applied to compute additional left/right
//         fluxes. Pointers to these functions are specified upon construction
//         of the solver and are used for non-standard fluxes (e.g. cosmic ray
//         energy/fluxes)
//
//       d This step is implemented in a virtual function implemented by a
//         subclass. The wavespeeds at the current interface is computed.
//         Then for each conserved (non-passive scalar) field, the array of
//         Riemann Flux (1C) at the current interface, is set equal to the
//         flux computed from the left/right reconstructed conserved values
//         (2B) and left/right fluxes (2C). This is achieved by iterating over
//         the indices of each array simultaneously.
//
//   - Use of the following field_lut struct:
//
//       - The calculation of fluxes and wave speeds requires random access of
//         specific fields. We also need to be able to iterate over the entries
//         of multiple arrays simultaneously (e.g. to accomplish part d, above)
//         Unlike Enzo, we wanted to avoid statically declaring which indices
//         correspond to which fields (adding additional sets of fields, like
//         internal energy and cosmic ray energy/fluxes becomes harder)
//
//       - We settled on using the field_lut struct as a lookup table. The
//         struct has members named for every quantity listed in FIELD_TABLE
//           - For a SCALAR, the member name directly matches the name in
//             column 1
//           - For a VECTOR, there are 3 members: {name}_i, {name}_j, {name}_k
//             ({name} cooresponds to the name appearing in column 1)
//         Each struct contains members named for all quantities in the table
//         (it includes conserved AND primitive quantites).
//
//       - Example: If we have an array of reconstructed primitives, wl, and
//         an instance of field_lut, prim_lut, that stores the indices of
//         primitives, then wl[prim_lut.density] and wl[prim_lut.pressure]
//         indicates the entries reserved for density and pressure (an
//         instance of field_lut storing indices for conserved quantities
//         would also indicate the index where density - since density is BOTH
//         conserved AND primitive)
//
//       - Given an instance of EnzoFieldConditions, the prepare_conserved_lut
//         function yields an initialized instance of field_lut and the
//         the length necessary for an array to hold values related to
//         conserved quantities. All members of field_lut corresponding to
//         conserved quantities are set equal to consectuive integers starting
//         from 0 (All members that don't correspond to conserved quantites
//         are set to -1).
//
//       - prepare_primitive_lut is analogous to prepare_conserved_lut except
//         it applies to primitive quantities
//     
//       - load_array_of_fields constructs an array of instances of EFlt3DArray
//         that correspond to reconstructed primitive fields, reconstructed
//         conserved fields OR flux fields. The function requires a pointer to
//         an instance of Block, the relevant initialized field_lut, the
//         length of the output field, an instance of grouping (group names
//         must match the relevant quantity names) and the direction along
//         which we are computing fluxes


// To check that the static methods of ImplStruct are appropriately defined,
// we follow https://stackoverflow.com/a/23133787
#define DEFINE_HAS_SIGNATURE(traitsName, funcName, signature)               \
  template <typename U>                                                     \
  class traitsName                                                          \
  {                                                                         \
  private:                                                                  \
    template<typename T, T> struct helper;                                  \
    template<typename T>                                                    \
    static std::uint8_t check(helper<signature, &funcName>*);               \
    template<typename T> static std::uint16_t check(...);                   \
  public:                                                                   \
    static                                                                  \
    constexpr bool value = sizeof(check<U>(0)) == sizeof(std::uint8_t);	    \
  }

DEFINE_HAS_SIGNATURE(has_calc_riemann_fluxes, T::calc_riemann_fluxes,
		     void (*)(const enzo_float[], const enzo_float[],
			      const enzo_float[], const enzo_float[],
			      const enzo_float[], const enzo_float[],
			      const field_lut, const field_lut, const int,
			      const bool, const enzo_float, const enzo_float,
			      const int, const int, const int, EFlt3DArray[],
			      enzo_float[]));
DEFINE_HAS_SIGNATURE(has_scratch_space_length, T::scratch_space_length,
		     int (*)(const int));

template <class ImplStruct>
class EnzoRiemannImpl : public EnzoRiemann
{
  /// @class    EnzoRiemannImpl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides implementation of approximate Riemann
  /// Solvers

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
  EnzoRiemannImpl(EnzoFieldConditions cond,
		  std::vector<std::string> &extra_passive_groups,
		  FluxFunctor** flux_funcs, int n_funcs);

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
	      Grouping &flux_group, Grouping &consl_group,
	      Grouping &consr_group, int dim, EnzoEquationOfState *eos,
	      int stale_depth);

protected : //methods
  
  /// Computes the fluxes for the passively advected quantites.
  void solve_passive_advection_(Block* block, Grouping &priml_group,
				Grouping &primr_group, Grouping &flux_group,
				EFlt3DArray &density_flux, int dim);

  /// computes fluxes for the basic mhd conserved quantities - density,
  /// momentum, energy, magnetic fields  
  void basic_mhd_fluxes_(const enzo_float prim[], const enzo_float cons[],
			 enzo_float fluxes[], const field_lut prim_lut,
			 const field_lut cons_lut);

protected: //attributes

  /// Conditions of the calculation used to initialized instances of field_lut
  EnzoFieldConditions conditions_;

  /// struct lookup-table that maps conserved field names to indices
  field_lut cons_lut_;

  /// struct lookup-table that maps primitive field names to indices
  field_lut prim_lut_;

  /// Number of conserved keys (fields) in cons_lut_
  int n_cons_keys_;

  /// Number of primitive keys (fields) in prim_lut_
  int n_prim_keys_;

  /// Names of the passively advected groups of fields (e.g. colors)
  std::vector<std::string> passive_groups_;

  /// number of flux functors
  int n_funcs_;

  /// array of pointers to functors used to compute fluxes
  FluxFunctor** flux_funcs_;

};

//----------------------------------------------------------------------

template <class ImplStruct>
EnzoRiemannImpl<ImplStruct>::EnzoRiemannImpl
(EnzoFieldConditions cond, std::vector<std::string> &passive_groups,
 FluxFunctor** flux_funcs, int n_funcs)
  : EnzoRiemann()
{
  conditions_ = cond;
  EnzoCenteredFieldRegistry registry;
  cons_lut_ = registry.prepare_cons_lut(conditions_, &n_cons_keys_);
  prim_lut_ = registry.prepare_prim_lut(conditions_, &n_prim_keys_);

  // construct vector of all groups of passively advected fields
  std::vector<std::string> passive_groups_ = passive_groups;

  n_funcs_ = n_funcs;
  flux_funcs_ = flux_funcs;
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
    // avoiding PUPing field_luts
    EnzoCenteredFieldRegistry registry;
    cons_lut_ = registry.prepare_cons_lut(conditions_, &n_cons_keys_);
    prim_lut_ = registry.prepare_prim_lut(conditions_, &n_prim_keys_);
  }

  p|passive_groups_;

  p|n_funcs_;

  if (p.isUnpacking()){
    flux_funcs_ = new FluxFunctor*[n_funcs_];
  }
  for (int i = 0; i<n_funcs_;i++){
    p|flux_funcs_[i];
  }
}

//----------------------------------------------------------------------

template <class ImplStruct>
void EnzoRiemannImpl<ImplStruct>::solve
(Block *block, Grouping &priml_group, Grouping &primr_group,
 Grouping &flux_group, Grouping &consl_group, Grouping &consr_group,
 int dim, EnzoEquationOfState *eos, int stale_depth)
{
  
  const bool barotropic = eos->is_barotropic();
  // if not barotropic then the following doesn't have to be reasonable
  const enzo_float isothermal_cs = eos->get_isothermal_sound_speed();
  // If barotropic, then the following doesn't have to be a reasonable value
  // (this will have to be tweaked when we introduce species that modify gamma)
  const enzo_float gamma = eos->get_gamma();

  // When barotropic equations of state are eventually introduced, all eos
  // dependencies should be moved up here
  ASSERT("EnzoRiemannImpl", "currently no support for barotropic eos",
	 !barotropic);
  
  EnzoCenteredFieldRegistry registry;
  EFlt3DArray *wl_arrays, *wr_arrays, *ul_arrays, *ur_arrays, *flux_arrays;

  wl_arrays = registry.load_array_of_fields(block, prim_lut_, n_prim_keys_,
					    priml_group, dim);
  wr_arrays = registry.load_array_of_fields(block, prim_lut_, n_prim_keys_,
					    primr_group, dim);
  ul_arrays = registry.load_array_of_fields(block, cons_lut_, n_cons_keys_,
					    consl_group, dim);
  ur_arrays = registry.load_array_of_fields(block, cons_lut_, n_cons_keys_,
					    consr_group, dim);
  flux_arrays = registry.load_array_of_fields(block, cons_lut_, n_cons_keys_,
					      flux_group, dim);

  enzo_float *wl, *wr, *Ul, *Ur, *Fl, *Fr;
  wl = new enzo_float[n_prim_keys_];  wr = new enzo_float[n_prim_keys_];
  Ul = new enzo_float[n_cons_keys_];  Ur = new enzo_float[n_cons_keys_];
  Fl = new enzo_float[n_cons_keys_];  Fr = new enzo_float[n_cons_keys_];

  enzo_float *scratch_space = NULL;
  int scratch_space_length = ImplStruct::scratch_space_length(n_cons_keys_);
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
	for (int field_ind=0; field_ind<n_prim_keys_; field_ind++){
	  wl[field_ind] = wl_arrays[field_ind](iz,iy,ix);
	  wr[field_ind] = wr_arrays[field_ind](iz,iy,ix);
	}
	for (int field_ind=0; field_ind<n_cons_keys_; field_ind++){
	  Ul[field_ind] = ul_arrays[field_ind](iz,iy,ix);
	  Ur[field_ind] = ur_arrays[field_ind](iz,iy,ix);
	}

	// compute the interface fluxes
	basic_mhd_fluxes_(wl, Ul, Fl, prim_lut_, cons_lut_);
	basic_mhd_fluxes_(wr, Ur, Fr, prim_lut_, cons_lut_);

	// iterate over the functors
	for (int i = 0; i<n_funcs_; i++){
	  (*(flux_funcs_[i]))(wl, Ul, Fl, prim_lut_, cons_lut_);
	  (*(flux_funcs_[i]))(wr, Ur, Fr, prim_lut_, cons_lut_);
	}


	// Now compute the Riemann Fluxes
	ImplStruct::calc_riemann_fluxes(Fl, Fr, wl, wr, Ul, Ur,
					prim_lut_, cons_lut_, n_cons_keys_,
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
void EnzoRiemannImpl<ImplStruct>::basic_mhd_fluxes_
(const enzo_float prim[], const enzo_float cons[], enzo_float fluxes[],
 const field_lut prim_lut, const field_lut cons_lut)
{
  // This assumes that MHD is included
  // This may be better handled by the EquationOfState
  enzo_float vi, vj, vk, p, Bi, Bj, Bk, etot, mag_pressure;
  vi = prim[prim_lut.velocity_i];
  vj = prim[prim_lut.velocity_j];
  vk = prim[prim_lut.velocity_k];
  p  = prim[prim_lut.pressure];
  Bi = prim[prim_lut.bfield_i];
  Bj = prim[prim_lut.bfield_j];
  Bk = prim[prim_lut.bfield_k];

  mag_pressure = mag_pressure_(prim, prim_lut);

  // Compute Fluxes
  enzo_float pi = cons[cons_lut.momentum_i];
  fluxes[cons_lut.density] = pi;

  // Fluxes for Mx, My, Mz
  fluxes[cons_lut.momentum_i] = pi*vi - Bi*Bi + p + mag_pressure;
  fluxes[cons_lut.momentum_j] = pi*vj - Bj*Bi;
  fluxes[cons_lut.momentum_k] = pi*vk - Bk*Bi;

  // Flux for etot
  etot = cons[cons_lut.total_energy];
  fluxes[cons_lut.total_energy] = ((etot + p + mag_pressure)*vi
				  - (Bi*vi + Bj*vj + Bk*vk)*Bi);

  // Fluxes for Bi,Bj,Bk
  fluxes[cons_lut.bfield_i] = 0;
  fluxes[cons_lut.bfield_j] = Bj*vi - Bi*vj;
  fluxes[cons_lut.bfield_k] = Bk*vi - Bi*vk;
}

#endif /* ENZO_ENZO_RIEMANN_IMPL_HPP */
