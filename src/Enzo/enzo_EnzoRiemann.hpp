
#ifndef ENZO_ENZO_RIEMANN_HPP
#define ENZO_ENZO_RIEMANN_HPP

#include <pup_stl.h>

// This class factors out the repeated code between different approximate
// Riemann Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF & Roe solvers).
// Between each of these solvers, only the code to compute the wave speeds and
// actual Riemann Fluxes changes.
// - if the virtual function calls are too slow, this implementation can be
//   refactored to use templates: EnzoRiemannHLLBase<S>. In that case, S would
//   refer to a helper class with a wave_speeds_ and riemann_flux_ helper
//   function. We could then name template specializations
//   (e.g. typedef EnzoRiemannHLLBase<HLLEHelper> EnzoRiemannHLLE)
// The code internally uses a struct as a lookup to map array indices to
// hydrodynamic quantities.
//
// To allow for easily adding additional fields with non-trivial flux
// calculations (e.g. cosmic rays or possibly internal energy), the
// constructor will accept an array of functors, subclassed from FluxFunctor
// - Again, if the functors are too slow, we could specify the functors as
//   variadic template arguments (allowing operator() to be inlined within the
//   for loop)
//
// To add an additional fields, add the field to FLUX_TABLE and modify the
// factory method
//
// Current state:
//  - support for passive scalars is not yet implemented (though there is a
//    spot carved out for it)


// Derivatives of FIELD_TABLE used by EnzoRiemann
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

class EnzoRiemann : public PUP::able
{
  /// @class    EnzoRiemann
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate approximate Riemann Solvers

public: // interface

  // Factory method for constructing the EnzoRiemann object
  // The signature of this method must be modified as additional physics get's
  // added
  static EnzoRiemann* construct_riemann(std::string solver);

  EnzoRiemann(EnzoFieldConditions cond,
	      std::vector<std::string> &extra_passive_groups,
	      FluxFunctor** flux_funcs, int n_funcs);

  /// Virtual destructor
  virtual ~EnzoRiemann();

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoRiemann);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemann (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void solve (Block *block, Grouping &priml_group, Grouping &primr_group,
	      Grouping &flux_group, Grouping &consl_group,
	      Grouping &consr_group, int dim, EnzoEquationOfState *eos);

protected : //methods

  // Compute the Riemann Fluxes (and wave_speeds)
  virtual void calc_riemann_fluxes_(const enzo_float flux_l[],
				    const enzo_float flux_r[],
				    const enzo_float prim_l[],
				    const enzo_float prim_r[],
				    const enzo_float cons_l[],
				    const enzo_float cons_r[],
				    const field_lut prim_lut,
				    const field_lut cons_lut,
				    const int n_keys,
				    EnzoEquationOfState *eos,
				    const int iz, const int iy, const int ix,
				    EFlt3DArray flux_arrays[]) =0;

  void solve_passive_advection_(Block* block,
				Grouping &priml_group, Grouping &primr_group,
				EFlt3DArray &density_flux, int dim)
  { /* This needs to be implemented */ }

  // computes the fast magnetosonic speed along dimension i
  enzo_float fast_magnetosonic_speed_(const enzo_float prim_vals[],
				      const field_lut prim_lut,
				      EnzoEquationOfState *eos);

  enzo_float mag_pressure_(const enzo_float prim_vals[],
			   const field_lut prim_lut);

  enzo_float sound_speed_(const enzo_float prim_vals[],
			  const field_lut prim_lut,
			  EnzoEquationOfState *eos)
  {
    return std::sqrt(eos->get_gamma()*prim_vals[prim_lut.pressure]/
		     prim_vals[prim_lut.density]);
  }

  void basic_mhd_fluxes_(const enzo_float prim[], const enzo_float cons[],
			 enzo_float fluxes[], const field_lut prim_lut,
			 const field_lut cons_lut);

protected: //attributes

  EnzoFieldConditions conditions_;
  field_lut cons_lut_;
  field_lut prim_lut_;
  int n_cons_keys_;
  int n_prim_keys_;

  std::vector<std::string> passive_groups_;
  // number of flux functors
  int n_funcs_;
  // array of pointers to functors used to compute fluxes
  FluxFunctor** flux_funcs_;
  enzo_float Dedner_Ch_;
};

#endif /* ENZO_ENZO_RIEMANN_HPP */
