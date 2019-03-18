#include <string>
#include <algorithm>
#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoRiemann* EnzoRiemann::construct_riemann(std::string solver)
{
  

  // In the future, need to construct vectors of extra group names and allocate
  // array of flux functors here

  // determine the type of solver to construct:
  
  // convert string to lower case (https://stackoverflow.com/a/313990)
  std::string formatted;
  std::transform(solver.begin(), solver.end(), formatted.begin(), std::tolower);

  ASSERT("EnzoRiemann", "The only allowed Solver is currently HLLE",
	 formatted == std::string("hlle"));

  EnzoRiemann* out = new EnzoRiemannHLLE();

  return out;
}

//----------------------------------------------------------------------

EnzoRiemann::EnzoRiemann(std::vector<std::string> &extra_scalar_groups,
			 std::vector<std::string> &extra_vector_groups,
			 std::vector<std::string> &extra_passive_groups,
			 FluxFunctor** flux_funcs, int n_funcs)
{
  // construct vectors of all cons/prim groups
  std::vector<std::string> d_cons_groups{"density", "momentum",
      "total_energy", "bfield"};
  std::vector<std::string> d_prim_groups{"density", "velocity",
      "pressure", "bfield"};
  combine_groups_(d_cons_groups, extra_scalar_groups, extra_vector_groups,
		  cons_groups_);
  combine_groups_(d_prim_groups, extra_scalar_groups, extra_vector_groups,
		  prim_groups_);

  // construct of all groups of passively advected fields
  std::vector<std::string> d_passive_groups{"species", "colors"};
  passive_groups_.reserve(d_passive_groups.size() +
			  extra_passive_groups.size());
  passive_groups_.insert(passive_groups_.end(), d_passive_groups.begin(),
			 d_passive_groups.end());
  passive_groups_.insert(passive_groups_.end(), extra_passive_groups.begin(),
			 extra_passive_groups.end());

  // We don't include passive_advected_groups in the main calculation so the
  // the values are not actually read into the maps
  n_keys_ = (8 + extra_scalar_groups.size() +
	     3*extra_scalar_groups.size());

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

  p|cons_groups_;
  p|prim_groups_;
  p|n_keys_;
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
  flt_map wl, wr, Ul, Ur, Fl, Fr;
  wl.reserve(n_keys_); wr.reserve(n_keys_);
  Ul.reserve(n_keys_); Ur.reserve(n_keys_);
  Fl.reserve(n_keys_); Fr.reserve(n_keys_);

  // Load in the fields
  array_map wl_arrays, wr_arrays, ul_arrays, ur_arrays, flux_arrays;
  wl_arrays.reserve(n_keys_); wr_arrays.reserve(n_keys_);
  ul_arrays.reserve(n_keys_); ur_arrays.reserve(n_keys_);
  flux_arrays.reserve(n_keys_);

  std::vector<std::string> prim_keys, cons_keys;
  load_fluid_fields_(block, wl_arrays, priml_group, prim_groups_, dim,
		     &prim_keys);
  load_fluid_fields_(block, wr_arrays, primr_group, prim_groups_, dim, NULL);
  load_fluid_fields_(block, ul_arrays, consl_group, cons_groups_, dim,
		     &cons_keys);
  load_fluid_fields_(block, ur_arrays, consr_group, cons_groups_, dim, NULL);
  load_fluid_fields_(block, flux_arrays, flux_group, cons_groups_, dim, NULL);

  // For Nearest-Neighbor, we care about interfaces starting at i+1/2
  // For PLM we only care about interfaces starting at i+3/2.
  //    (left interface value at i+1/2 is set to 0)
  // For consistency, always start at i+1/2
  // Iteration limits are generalized for 2D and 3D arrays

  for (int iz=0; iz<flux_arrays["density"].shape(0); iz++) {
    for (int iy=0; iy<flux_arrays["density"].shape(1); iy++) {
      for (int ix=0; ix<flux_arrays["density"].shape(2); ix++) {

	// get the fluid fields
	for (unsigned int field_ind=0; field_ind<n_keys_; field_ind++){
	  std::string key = prim_keys[field_ind];
	  wl[key] = wl_arrays[key](iz,iy,ix);
	  wr[key] = wr_arrays[key](iz,iy,ix);

	  key = cons_keys[field_ind];
	  Ul[key] = ul_arrays[key](iz,iy,ix);
	  Ur[key] = ur_arrays[key](iz,iy,ix);
	}

	// compute the interface fluxes
	basic_mhd_fluxes_(wl, Ul, Fl);
	basic_mhd_fluxes_(wr, Ur, Fr);

	// iterate over the functors
	for (int i = 0; i<n_funcs_; i++){
	  (*(flux_funcs_[i]))(wl, Ul, Fl);
	  (*(flux_funcs_[i]))(wr, Ur, Fr);
	}

	// Now compute the Riemann Fluxes
	calc_riemann_fluxes_(Fl, Fr, wl, wr, Ul, Ur, cons_keys,
			     n_keys, eos, iz, iy, ix, flux_arrays);

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
	  flux_arrays["phi"] *= (Dedner_Ch_*Dedner_Ch_);
	}
      }
    }
  }
  solve_passive_advection_(block, wl, wr, flux_arrays["density"], dim);
}

//======================================================================

void EnzoRiemann::combine_groups_(std::vector<std::string> &default_g,
				  std::vector<std::string> &extra_scalar_g,
				  std::vector<std::string> &extra_vector_g,
				  std::vector<std::string> &combined_g)
{
  int n_groups = 4 + extra_scalar_groups.size() + extra_vector_groups.size();
  combined_g.reserve(n_groups);
  combined_g.insert(combined_g.end(), default_g.begin(), default_g.end());
  combined_g.insert(combined_g.end(), extra_scalar_g.begin(),
		    extra_scalar_g.end());
  combined_g.insert(combined_g.end(), extra_vector_g.begin(),
		    extra_vector_g.end());
}

//----------------------------------------------------------------------

void EnzoRiemann::load_fluid_fields_(Block *block, array_map &arrays,
				     Grouping &grouping,
				     std::vector<std::string> &group_names,
				     int dim,
				     std::vector<std::string> *key_names)
{
  // If key_names != NULL, insert the name of each key
  EnzoPermutedCoordinates coord(dim);
  EnzoFieldArrayFactory array_factory(block);
 
  for (std::size_t i; i<group_names.size(); i++){
    std::string group_name = group_names[i];
    int group_size = grouping.size(group_name);
    ASSERT("EnzoRiemannHLLBase",
	   "Implementation requires non-passive groups to have 1 or 3 fields",
	   group_size == 1 || group_size == 3);

    for (int num = 0; j<group_size; j++){
      int field_ind;
      std::string key;
      if (group_size == 1){
	key = group_name;
	field_ind = 0;
      } else if (num == 0) {
	key = group_name + std::string("_i");
	field_ind = coord.i_axis();
      } else if (num == 1) {
	key = group_name + std::string("_j");
	field_ind = coord.j_axis();
      } else {
	key = group_name + std::string("_k");
	field_ind = coord.k_axis();
      }
      arrays[key] = array_factory.reconstructed_field(grouping, group_name,
						      field_ind, dim);
      if (key_names != NULL){
	key_names.push_back(key);
      }
    }      
  }
}

//----------------------------------------------------------------------

enzo_float EnzoRiemann::fast_magnetosonic_speed_(const flt_map &prim_vals)
{
  enzo_float cs2 = std::pow(sound_speed(prim_vals),2);
  enzo_float B2 = (prim_vals["bfield_i"]*prim_vals["bfield_i"] +
		   prim_vals["bfield_j"]*prim_vals["bfield_j"] +
		   prim_vals["bfield_k"]*prim_vals["bfield_k"]);
  enzo_float va2 = B2/prim_vals["density"];
  enzo_float cos2 = prim_vals["bfield_i"]*prim_vals["bfield_i"] / B2;
  return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
					  4.*cs2*va2*cos2)));
}

//----------------------------------------------------------------------

enzo_float EnzoRiemann::mag_pressure_(const flt_map &prim_vals)
{
  return 0.5 * (prim_vals["bfield_i"]*prim_vals["bfield_i"]+
		prim_vals["bfield_j"]*prim_vals["bfield_j"]+
		prim_vals["bfield_k"]*prim_vals["bfield_k"]);
}

//----------------------------------------------------------------------

void EnzoRiemann::basic_mhd_fluxes_(const flt_map &prim, const flt_map &cons,
				    flt_map &fluxes)
{
  // This assumes that MHD is included
  // This may be better handled by the EquationOfState
  enzo_float vi, vj, vk, p, Bi, Bj, Bk, etot, mag_pressure;
  vi = prim["velocity_i"];
  vj = prim["velocity_j"];
  vk = prim["velocity_k"];
  p  = prim["pressure"];
  Bi = prim["bfield_i"];
  Bj = prim["bfield_j"];
  Bk = prim["bfield_k"];

  mag_pressure = mag_pressure_(prim);

  // Compute Fluxes
  fluxes["density"] = cons["momentum_i"];

  // Fluxes for Mx, My, Mz
  fluxes["momentum_i"] = fluxes["density"]*vi - Bi*Bi + p + mag_pressure;
  fluxes["momentum_j"] = fluxes["density"]*vj - Bj*Bi;
  fluxes["momentum_k"] = fluxes["density"]*vk - Bk*Bi;

  // Flux for etot
  etot = cons["total_energy"];
  fluxes["total_energy"] = ((etot + p + mag_pressure)*vi
			    - (Bi*vi + Bj*vj + Bk*vk)*Bi);

  // Fluxes for Bi,Bj,Bk
  fluxes["bfield_i"] = 0;
  fluxes["bfield_j"] = Bj*vi - Bi*vj;
  fluxes["bfield_k"] = Bk*vi - Bi*vk;
}
