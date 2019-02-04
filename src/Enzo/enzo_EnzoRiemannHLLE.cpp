#include "cello.hpp"
#include "enzo.hpp"

// Some Notes
// Stone+ (08) has a typo in eqn 52: (q_i-q_{i-1}) should be (q_R-q_L)
//    The modified formula is actually used in Athena++

// Implementation Notes:
//   - Based off of the HLLE Riemann Solver described in Stone+08
//   - Currently only supports adiabatic fluids
//        To fix this, may want to shift the responsibility of computing
//        eigenvalues of the Roe matrix to EnzoEquationOfState
//        Benefit: Decouples Riemann Solver from fluid type
//                 Probably could reuse code in a Roe Solver
//   - This class is primarily implemented to allow for easy extension to
//     incorporate CRs (just need to swap out the wave_speeds_ helper method)
//  - If the need arises, have to modify to handle dual energy formalism

typedef std::unordered_map<std::string,EFlt3DArray> array_map;

void add_entry_(EnzoFieldArrayFactory &af, array_map &arrays,
		std::string key, Grouping &grouping, std::string group_name,
		int index, int dim)
{
  arrays[key] = af.load_temp_interface_grouping_field(grouping, group_name,
						      index, dim != 0,
						      dim != 1, dim != 2);
}

void load_fluid_fields_(Block *block, array_map &arrays, Grouping &grouping,
			bool primitive, int dim)
{
  int i = dim;
  int j = (dim+1)%3;
  int k = (dim+2)%3;

  EnzoFieldArrayFactory af(block);
  
  // Load density
  add_entry_(af, arrays, "density", grouping, "density", 0, dim);

  // load velocity/momentum and pressure/total energy
  if (primitive){
    add_entry_(af, arrays, "velocity_i", grouping, "velocity", i, dim);
    add_entry_(af, arrays, "velocity_j", grouping, "velocity", j, dim);
    add_entry_(af, arrays, "velocity_k", grouping, "velocity", k, dim);
    add_entry_(af, arrays, "pressure", grouping, "pressure", 0, dim);
  } else {
    add_entry_(af, arrays, "momentum_i", grouping, "momentum", i, dim);
    add_entry_(af, arrays, "momentum_j", grouping, "momentum", j, dim);
    add_entry_(af, arrays, "momentum_k", grouping, "momentum", k, dim);
    add_entry_(af, arrays, "total_energy", grouping, "total_energy", 0, dim);
  }

  add_entry_(af, arrays, "bfield_i", grouping, "bfield", i, dim);
  add_entry_(af, arrays, "bfield_j", grouping, "bfield", j, dim);
  add_entry_(af, arrays, "bfield_k", grouping, "bfield", k, dim);
}



//----------------------------------------------------------------------

void EnzoRiemannHLLE::solve (Block *block, Grouping &priml_group,
			     Grouping &primr_group, Grouping &flux_group,
			     Grouping &consl_group, Grouping &consr_group,
			     int dim, EnzoEquationOfState *eos)
{

  // it would probably make sense for prim_keys and cons_keys to be instance
  // variables of the Riemann object
  std::vector<std::string> prim_keys{"density", "velocity_i", "velocity_j",
      "velocity_k", "pressure", "bfield_i", "bfield_j", "bfield_k"};
  std::vector<std::string> cons_keys{"density", "momentum_i", "momentum_j",
      "momentum_k", "total_energy", "bfield_i", "bfield_j", "bfield_k"};
  unsigned int length = 8;

  flt_map wl, wr, Ul, Ur, Fl, Fr;
  wl.reserve(length); wr.reserve(length);
  Ul.reserve(length); Ur.reserve(length);
  Fl.reserve(length); Fr.reserve(length);
  enzo_float bp, bm;
  
  //int nspecies = priml_groups.size("species");
  //int ncolors = priml_groups.size("colors");

  // Load in the fields
  array_map wl_arrays, wr_arrays, ul_arrays, ur_arrays, flux_arrays;
  wl_arrays.reserve(length); wr_arrays.reserve(length);
  ul_arrays.reserve(length); ur_arrays.reserve(length);
  flux_arrays.reserve(length);

  load_fluid_fields_(block, wl_arrays, priml_group, true, dim);
  load_fluid_fields_(block, wr_arrays, primr_group, true, dim);
  load_fluid_fields_(block, ul_arrays, consl_group, false, dim);
  load_fluid_fields_(block, ur_arrays, consr_group, false, dim);
  load_fluid_fields_(block, flux_arrays, flux_group, false, dim);

  // For Nearest-Neighbor, we care about interfaces starting at i+1/2
  // For PLM we only care about interfaces starting at i+3/2.
  //    (left interface value at i+1/2 is set to 0)
  // For consistency, always start at i+1/2
  // Iteration limits are generalized for 2D and 3D arrays

  for (int iz=0; iz<flux_arrays["density"].length_dim2(); iz++) {
    for (int iy=0; iy<flux_arrays["density"].length_dim1(); iy++) {
      for (int ix=0; ix<flux_arrays["density"].length_dim0(); ix++) {

	// get the fluid fields
	for (unsigned int field_ind=0; field_ind<length; field_ind++){
	  std::string key = prim_keys[field_ind];
	  wl[key] = wl_arrays[key](iz,iy,ix);
	  wr[key] = wr_arrays[key](iz,iy,ix);

	  key = cons_keys[field_ind];
	  Ul[key] = ul_arrays[key](iz,iy,ix);
	  Ur[key] = ur_arrays[key](iz,iy,ix);
	}

	// Compute magnetic pressure
	enzo_float mag_pressure_l, mag_pressure_r;
	mag_pressure_l = eos->mag_pressure_from_primitive(wl);
	mag_pressure_r = eos->mag_pressure_from_primitive(wr);

	// compute the wavespeeds
	wave_speeds_(wl, wr, Ul, Ur, mag_pressure_l, mag_pressure_r,
		     eos, &bp, &bm);

	// compute the interface fluxes
	interface_flux_(wl, Ul, Fl, mag_pressure_l);
	interface_flux_(wr, Ur, Fr, mag_pressure_r);

	// Now compute the Riemann Flux
	for (unsigned int field_ind=0; field_ind<length; field_ind++){
	  std::string key = cons_keys[field_ind];

	  flux_arrays[key](iz,iy,ix) =
	    ((bp*Fl[key] - bm*Fr[key]  +
	      (Ur[key] - Ul[key])*bp*bm)/ (bp - bm));

	  ASSERT("EnzoRiemannHLLE",
	  	 "There is a NaN or inf\n",
	  	 std::isfinite(flux_arrays[key](iz,iy,ix)));
	}

	// Deal with Species and colors
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoRiemannHLLE::wave_speeds_(flt_map &wl, flt_map &wr, flt_map &Ul,
				   flt_map &Ur, enzo_float mag_p_l,
				   enzo_float mag_p_r,
				   EnzoEquationOfState *eos,
				   enzo_float *bp, enzo_float *bm)
{
  // Calculate wavespeeds as specified by S4.3.1 of Stone+08
  // if LM and L0 are max/min eigenvalues of Roe's matrix:
  //       bp = max(LM, vi_r + cf_r, 0) (eqn 53)
  //       bm = max(L0, vi_l - cf_l, 0) (eqn 54)
  // vi_r and vi_l are velocities of left and right states along ith dimension
  // cf_r and cf_l are the fast magnetosonic speeds in the left and right states

  // First, compute left and right speeds
  //    left_speed = vi_l - cf_l;      right_speed = vi_r + cf_r
  enzo_float left_speed = wl["velocity_i"] - eos->fast_magnetosonic_speed(wl);
  enzo_float right_speed = wr["velocity_i"] + eos->fast_magnetosonic_speed(wr);

  // Next, compute min max eigenvalues
  //     - per eqn B17 these are Roe averaged velocity in the ith direction
  //       minus/plus Roe averaged fast magnetosonic wavespeed
  //     - Will probably offload this to a method of EnzoEquationOfState
  enzo_float sqrtrho_l = std::sqrt(wl["density"]);
  enzo_float sqrtrho_r = std::sqrt(wr["density"]);
  enzo_float coef = 1.0/(sqrtrho_l + sqrtrho_r);

  // density and velocity
  enzo_float rho_roe = sqrtrho_l*sqrtrho_r;
  enzo_float vi_roe = (sqrtrho_l * wl["velocity_i"] +
		       sqrtrho_r * wr["velocity_i"])*coef;
  enzo_float vj_roe = (sqrtrho_l * wl["velocity_j"] +
		       sqrtrho_r * wr["velocity_j"])*coef;
  enzo_float vk_roe = (sqrtrho_l * wl["velocity_k"] +
		       sqrtrho_r * wr["velocity_k"])*coef;

  // enthalpy:  h = (etot + thermal pressure + magnetic pressure) / rho
  enzo_float h_l = ((Ul["total_energy"] + wl["pressure"] + mag_p_l) /
		    wl["density"]);
  enzo_float h_r = ((Ur["total_energy"] + wr["pressure"] + mag_p_r) /
		    wr["density"]);
  enzo_float h_roe = (sqrtrho_l * h_l + sqrtrho_r * h_r) * coef;

  // Magnetic Fields (formulas are different pattern from above)
  enzo_float bi_roe = wl["bfield_i"]; // equal to wr["bfield_i"]
  enzo_float bj_roe = (sqrtrho_l * wr["bfield_j"] +
		       sqrtrho_r * wl["bfield_j"])*coef;
  enzo_float bk_roe = (sqrtrho_l * wr["bfield_k"] +
		       sqrtrho_r * wl["bfield_k"])*coef;


  // Calculate fast magnetosonic speed for Roe-averaged quantities (eqn B18)
  enzo_float gamma_prime = eos->get_gamma()-1.;
  enzo_float x_prime = ((std::pow(wl["bfield_j"]-wr["bfield_j"],2) +
			 std::pow(wl["bfield_k"]-wr["bfield_k"],2) )
			* 0.5 * (gamma_prime-1) * coef);
  enzo_float y_prime = ((gamma_prime-1)*(wl["density"]+wr["density"])
			*0.5/rho_roe);

  enzo_float v_roe2 = vi_roe*vi_roe + vj_roe*vj_roe + vk_roe*vk_roe;
  enzo_float b_roe2 = bi_roe*bi_roe + bj_roe*bj_roe + bk_roe*bk_roe;

  // Need tidle_a^2
  enzo_float tilde_a2 = (gamma_prime * (h_roe - 0.5* v_roe2 - b_roe2/rho_roe)
			 - x_prime);

  // Need Roe averaged square Alfven speed (along ith dim and total magnitude)
  enzo_float tilde_vai2 = bi_roe * bi_roe / rho_roe;
  enzo_float tilde_va2 = (tilde_vai2 + (gamma_prime - y_prime) *
			  (bj_roe * bj_roe + bk_roe * bk_roe) / rho_roe);
  //CkPrintf("roe_cs2 = %.15lf\n",tilde_a2);
  //CkPrintf("roe_va_i2 = %.15lf, roe_va2 = %.15lf\n",tilde_vai2, tilde_va2);

  enzo_float cfast = std::sqrt(0.5 * (tilde_a2 + tilde_va2 +
				      std::sqrt(std::pow(tilde_a2 + tilde_va2,
							 2) -
						4 * tilde_a2 * tilde_vai2)));

  //CkPrintf("vi_roe = %.15lf, cfast = %.15lf\n",vi_roe, cfast);
  //CkPrintf("left_speed = %.15lf, right_speed = %.15lf\n",
  //left_speed, right_speed);

  *bp = std::fmax(std::fmax(vi_roe + cfast, right_speed), 0.);
  *bm = std::fmin(std::fmin(vi_roe - cfast, left_speed),  0.);

  //CkPrintf("bp = %.15lf, bm = %.15lf\n",*bp,*bm);
  //fflush(stdout);
  //ASSERT("EnzoRiemannHLLE", "This is to stop the code!",false);
}



void EnzoRiemannHLLE::interface_flux_(flt_map &prim, flt_map &cons,
				      flt_map &fluxes, enzo_float mag_pressure)
{
  // can simplify things, like in Athena++ (described in Toro)
  // they use a modified interface flux form (it involves fewer steps)

  // This assumes that MHD is included
  // This may be better handled by the EquationOfState
  enzo_float vi, vj, vk, p, Bi, Bj, Bk, etot;
  vi = prim["velocity_i"];
  vj = prim["velocity_j"];
  vk = prim["velocity_k"];
  p  = prim["pressure"];
  Bi = prim["bfield_i"];
  Bj = prim["bfield_j"];
  Bk = prim["bfield_k"];

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
