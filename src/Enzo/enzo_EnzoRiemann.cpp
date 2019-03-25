#include <string>
#include <algorithm>
#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoRiemann* EnzoRiemann::construct_riemann(std::string solver)
{

  // Setup the EnzoFieldConditions struct
  EnzoFieldConditions cond;
  cond.hydro = true;
  cond.MHD = true;

  // Generate vector of group names of passive advected scalars
  std::vector<std::string> passive_groups{"species", "colors"};

  // In the future, allocate array of flux functors here
  FluxFunctor** flux_funcs = NULL;
  int n_funcs = 0;
    

  // determine the type of solver to construct:
  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted(solver.size(), ' ');
  std::transform(solver.begin(), solver.end(), formatted.begin(),
		 ::tolower);
  EnzoRiemann* out;

  if (formatted == std::string("hlle")){
    out = new EnzoRiemannHLLE(cond, passive_groups, flux_funcs, n_funcs);
  } else if (formatted == std::string("hlld")){
    out = new EnzoRiemannHLLD(cond, passive_groups, flux_funcs, n_funcs);
  } else {
    ASSERT("EnzoRiemann", "The only allowed solvers are HLLE & HLLD", false);
    out = NULL;  // Deals with compiler warning
  }

  return out;
}

//----------------------------------------------------------------------

EnzoRiemann::EnzoRiemann(EnzoFieldConditions cond,
			 std::vector<std::string> &passive_groups,
			 FluxFunctor** flux_funcs, int n_funcs)
{
  conditions_ = cond;
  cons_lut_ = prepare_conserved_lut(conditions_, &n_cons_keys_);
  prim_lut_ = prepare_primitive_lut(conditions_, &n_prim_keys_);
  
  // construct of all groups of passively advected fields
  std::vector<std::string> passive_groups_ = passive_groups;

  n_funcs_ = n_funcs;
  flux_funcs_ = flux_funcs;

  // For now we are not using Dedner Divergence cleaning
  Dedner_Ch_ = 0.;
}

//----------------------------------------------------------------------

EnzoRiemann::~EnzoRiemann()
{
  for (int i = 0; i < n_funcs_; i++){
    delete flux_funcs_[i];
  }
  delete[] flux_funcs_;
}

//----------------------------------------------------------------------

void EnzoRiemann::pup (PUP::er &p)
{
  PUP::able::pup(p);

  p|conditions_;
  if (p.isUnpacking()){
    // avoiding PUPing field_luts
    cons_lut_ = prepare_conserved_lut(conditions_, &n_cons_keys_);
    prim_lut_ = prepare_primitive_lut(conditions_, &n_prim_keys_);
  }

  p|passive_groups_;

  p|n_funcs_;

  if (p.isUnpacking()){
    flux_funcs_ = new FluxFunctor*[n_funcs_];
  }
  for (int i = 0; i<n_funcs_;i++){
    p|flux_funcs_[i];
  }

  p|Dedner_Ch_;
}

//----------------------------------------------------------------------

void EnzoRiemann::solve (Block *block, Grouping &priml_group,
			 Grouping &primr_group, Grouping &flux_group,
			 Grouping &consl_group, Grouping &consr_group, int dim,
			 EnzoEquationOfState *eos)
{
  // In the future, this stuff up here will be handled at construction

  EFlt3DArray *wl_arrays, *wr_arrays, *ul_arrays, *ur_arrays, *flux_arrays;


  wl_arrays = load_array_of_fields(block, prim_lut_, n_prim_keys_, priml_group,
				   dim);
  wr_arrays = load_array_of_fields(block, prim_lut_, n_prim_keys_, primr_group,
				   dim);
  ul_arrays = load_array_of_fields(block, cons_lut_, n_cons_keys_, consl_group,
				   dim);
  ur_arrays = load_array_of_fields(block, cons_lut_, n_cons_keys_, consr_group,
				   dim);
  flux_arrays = load_array_of_fields(block, cons_lut_, n_cons_keys_,
				     flux_group, dim);

  enzo_float *wl, *wr, *Ul, *Ur, *Fl, *Fr;
  wl = new enzo_float[n_prim_keys_];  wr = new enzo_float[n_prim_keys_];
  Ul = new enzo_float[n_cons_keys_];  Ur = new enzo_float[n_cons_keys_];
  Fl = new enzo_float[n_cons_keys_];  Fr = new enzo_float[n_cons_keys_];
  

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
	calc_riemann_fluxes_(Fl, Fr, wl, wr, Ul, Ur, prim_lut_, cons_lut_,
			     n_cons_keys_, eos, iz, iy, ix, flux_arrays);

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

  solve_passive_advection_(block, priml_group, primr_group,
			   flux_arrays[cons_lut_.density], dim);

  delete[] wl; delete[] wr;
  delete[] Ul; delete[] Ur;
  delete[] Fl; delete[] Fr;
  delete[] wl_arrays; delete[] wr_arrays;
  delete[] ul_arrays; delete[] ur_arrays;
  delete[] flux_arrays;
}

//======================================================================


enzo_float EnzoRiemann::fast_magnetosonic_speed_(const enzo_float prim_vals[],
						 const field_lut prim_lut,
						 EnzoEquationOfState *eos)
{
  enzo_float bi = prim_vals[prim_lut.bfield_i];
  enzo_float bj = prim_vals[prim_lut.bfield_j];
  enzo_float bk = prim_vals[prim_lut.bfield_k];
  
  enzo_float cs2 = std::pow(sound_speed_(prim_vals, prim_lut, eos),2);
  enzo_float B2 = (bi*bi + bj*bj + bk *bk);
  enzo_float va2 = B2/prim_vals[prim_lut.density];
  enzo_float cos2 = bi*bi / B2;
  return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
					  4.*cs2*va2*cos2)));
}

//----------------------------------------------------------------------

enzo_float EnzoRiemann::mag_pressure_(const enzo_float prim_vals[],
				      const field_lut prim_lut)
{
  enzo_float bi = prim_vals[prim_lut.bfield_i];
  enzo_float bj = prim_vals[prim_lut.bfield_j];
  enzo_float bk = prim_vals[prim_lut.bfield_k];
  return 0.5 * (bi*bi + bj*bj + bk *bk);
}

//----------------------------------------------------------------------

void EnzoRiemann::basic_mhd_fluxes_(const enzo_float prim[],
				    const enzo_float cons[],
				    enzo_float fluxes[],
				    const field_lut prim_lut,
				    const field_lut cons_lut)
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
