
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
//
// The code internally uses an unordered map to keep track of the hydrodynamic
// properties at a given cell in the mesh.
//typedef std::unordered_map<std::string,enzo_float> flt_map;
typedef std::unordered_map<std::string,EFlt3DArray> array_map;
// If the dynamic properties of a map slow the implementation too much, we
// could switch to an array based system (to determine ordering of values,
// templates can be used to calculate hash functions for strings at compile
// time)
//
// To allow for easily adding additional fields which require non-trivial
// flux calculations (e.g. cosmic rays or possibly internal energy),
// the constructor will accept an array of functors, subclassed from FluxFunctor
// Again, if the functors are too slow, we could specify the functors as
// variadic (allowing operator() to be inlined within the for loop)
//
// Current state:
//  - support for passive scalars is not yet implemented (though there is a
//    spot carved out for it)
//  -

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

  virtual void operator()(const flt_map &cons, const flt_map &prim,
			  flt_map &flux)=0;
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
					

  EnzoRiemann(std::vector<std::string> &extra_scalar_groups,
	      std::vector<std::string> &extra_vector_groups,
	      std::vector<std::string> &extra_passive_groups,
	      FluxFunctor** flux_funcs, int n_funcs)
  {
    initializer_(extra_scalar_groups, extra_vector_groups,
		 extra_passive_groups, flux_funcs, n_funcs);
  }

  EnzoRiemann()
  {
    // Empty vector
    std::vector<std::string> temp;
    initializer_(temp, temp, temp, NULL,0);
  };

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
  virtual void calc_riemann_fluxes_(flt_map &flux_l, flt_map &flux_r,
				    flt_map &prim_l, flt_map &prim_r,
				    flt_map &cons_l, flt_map &cons_r,
				    std::vector<std::string> &cons_keys,
				    std::size_t n_keys,
				    EnzoEquationOfState *eos,
				    const int iz, const int iy, const int ix,
				    array_map &flux_arrays) =0;

  void solve_passive_advection_(Block* block,
				Grouping &priml_group, Grouping &primr_group,
				EFlt3DArray &density_flux, int dim)
  { /* This needs to be implemented */ }

  void combine_groups_(std::vector<std::string> &default_g,
		       std::vector<std::string> &extra_scalar_g,
		       std::vector<std::string> &extra_vector_g,
		       std::vector<std::string> &combined_g);

  void load_fluid_fields_(Block *block, array_map &arrays, Grouping &grouping,
			  std::vector<std::string> &group_names, int dim,
			  std::vector<std::string> *key_names);

  // computes the fast magnetosonic speed along dimension i
  enzo_float fast_magnetosonic_speed_(const flt_map &prim_vals,
				      EnzoEquationOfState *eos);

  enzo_float mag_pressure_(const flt_map &prim_vals);
  
  enzo_float sound_speed_(const flt_map &prim_vals, EnzoEquationOfState *eos)
  {
    return std::sqrt(eos->get_gamma()*prim_vals.at("pressure")/
		     prim_vals.at("density"));
  }

  void basic_mhd_fluxes_(const flt_map &prim, const flt_map &cons,
			 flt_map &fluxes);
private:
  void initializer_(std::vector<std::string> &extra_scalar_groups,
		    std::vector<std::string> &extra_vector_groups,
		    std::vector<std::string> &extra_passive_groups,
		    FluxFunctor** flux_funcs, int n_funcs);
    
protected: //attributes
  std::vector<std::string> cons_groups_;
  std::vector<std::string> prim_groups_;
  std::size_t n_keys_;
  std::vector<std::string> passive_groups_;
  // number of flux functors
  int n_funcs_;
  // array of pointers to functors used to compute fluxes
  FluxFunctor** flux_funcs_;
  enzo_float Dedner_Ch_;
};

#endif /* ENZO_ENZO_RIEMANN_HPP */
