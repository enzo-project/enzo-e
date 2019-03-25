
#ifndef ENZO_ENZO_RIEMANN_HPP
#define ENZO_ENZO_RIEMANN_HPP

#include <pup_stl.h>
#include "enzo_EnzoRiemannFields.hpp"

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
