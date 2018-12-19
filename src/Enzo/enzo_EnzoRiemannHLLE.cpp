#include "cello.hpp"
#include "enzo.hpp"

// Some Notes
// Stone+ (08) has a typo in eqn 52: (q_i-q_{i-q}) should be (q_R-q_L)
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
//   - Unless otherwise noted, the indices in wl, wr, Fl, and Fr are
//     associated with quantities as follows (values in parenthesis indicate
//     the conserved quantity associated with flux values):
//       0 <-> density
//       1 <-> velocity_i (momentum_i)
//       2 <-> velocity_j (momentum_j)
//       3 <-> velocity_k (momentum_k)
//       4 <-> pressure   (fluid total energy density)
//       5 <-> b-field_i
//       6 <-> b-field_j
//       7 <-> b-field_k
//    Subclasses will probably do the same but add in extra fields at the end
//  - If the need arises, have to modify to handle dual energy formalism
//  - Going to replace the above with some kind of map (probably an unordered
//    map?)

void load_fluid_fields_(Field *field, std::vector<enzo_float*> *field_arrays,
			Grouping *grouping, bool primitive, int dim)
{
  int i = dim;
  int j = (dim+1)%3;
  int k = (dim+2)%3;

  // Load density
  field_arrays->push_back(load_grouping_field_(field, grouping, "density", 0));

  // load velocity/momentum
  std::string group_name = (primitive) ? "velocity" : "momentum";
  field_arrays->push_back(load_grouping_field_(field, grouping, group_name, i));
  field_arrays->push_back(load_grouping_field_(field, grouping, group_name, j));
  field_arrays->push_back(load_grouping_field_(field, grouping, group_name, k));

  // load pressure/total_energy
  group_name = (primitive) ? "pressure" : "total_energy";
  field_arrays->push_back(load_grouping_field_(field, grouping, group_name, 0));

  // load bfields
  field_arrays->push_back(load_grouping_field_(field, grouping, "bfield", i));
  field_arrays->push_back(load_grouping_field_(field, grouping, "bfield", j));
  field_arrays->push_back(load_grouping_field_(field, grouping, "bfield", k));
}


//----------------------------------------------------------------------

void EnzoRiemannHLLE::solve (Block *block, Grouping &priml_group,
			     Grouping &primr_group, Grouping &flux_group,
			     int dim, EnzoEquationOfState *eos)
{

  // it would probably make sense for length to be an instance variable of
  // Riemann object
  const int length = 8;
  // The following may not be legal
  enzo_float wl[length], wr[length], Ul[length], Ur[length];
  enzo_float Fl[length], Fr[length];
  enzo_float bp, bm;

  //int nspecies = priml_groups.size("species");
  //int ncolors = priml_groups.size("colors");

  // Declare the block and get field object
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // Load in the fields
  std::vector<enzo_float*> wl_arrays;
  std::vector<enzo_float*> wr_arrays;
  std::vector<enzo_float*> flux_arrays;

  load_fluid_fields_(&field, &wl_arrays, &priml_group, true, dim);
  load_fluid_fields_(&field, &wr_arrays, &primr_group, true, dim);
  load_fluid_fields_(&field, &flux_arrays, &flux_group, false, dim);

  // get integration limits 
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];

  // these values are for cell-centered grids. The need to be modified so that
  // along dim, the integration limits are face-centered
  switch(dim){
  case(0): mx++;
  case(1): my++;
  case(2): mz++;
  }

  // For PLM we only care about fluxes for the third cell in
  // For Nearest-Neighbor, we care about the second cell in
  // since mx, my, and mz are face-centered, the following should be correct
  for (int iz=1; iz<mz-1; iz++) {
    for (int iy=1; iy<my-1; iy++) {
      for (int ix=1; ix<mx-1; ix++) {
	// compute the index
	int i = ix + mx*(iy + my*iz);

	// get the primitive fields
	for (int field_ind=0; field_ind<length; field_ind++){
	  wl[field_ind] = wl_arrays[field_ind][i];
	  wr[field_ind] = wr_arrays[field_ind][i];
	}

	// compute the conservatives
	eos->conservative_from_primitive(wl, Ul);
	eos->conservative_from_primitive(wr, Ur);

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
	for (int field_ind=0; field_ind<length; field_ind++){
	  // fill in the following line properly
	  flux_arrays[field_ind][i] = (((bp*Fl[field_ind] - bm*Fr[field_ind])
				        / (bp - bm)) +
				       ((Ul[field_ind] - Ur[field_ind])*bp*bm
					/(bp - bm)));
	}

	// Deal with Species and colors
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoRiemannHLLE::wave_speeds_ (enzo_float *wl, enzo_float *wr,
				    enzo_float *Ul, enzo_float *Ur,
				    enzo_float mag_p_l, enzo_float mag_p_r,
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
  enzo_float left_speed = wl[1] - eos->fast_magnetosonic_speed(wl);
  enzo_float right_speed = wr[1] + eos->fast_magnetosonic_speed(wr);


  // Next, compute min max eigenvalues
  //     - per eqn B17 these are Roe averaged velocity in the ith direction
  //       minus/plus Roe averaged fast magnetosonic wavespeed
  //     - Will probably offload this to a method of EnzoEquationOfState
  enzo_float sqrtrho_l = std::sqrt(wl[0]);
  enzo_float sqrtrho_r = std::sqrt(wr[0]);
  enzo_float coef = 1.0/(sqrtrho_l + sqrtrho_r);

  // density and velocity
  enzo_float rho_roe = sqrtrho_l*sqrtrho_r;
  enzo_float vi_roe = (sqrtrho_l * wl[1] + sqrtrho_r * wr[1])*coef;
  enzo_float vj_roe = (sqrtrho_l * wl[2] + sqrtrho_r * wr[2])*coef;
  enzo_float vk_roe = (sqrtrho_l * wl[3] + sqrtrho_r * wr[3])*coef;

  // enthalpy:  h = (etot + thermal pressure + magnetic pressure) / rho
  enzo_float h_l = (Ul[4] + wl[4] + mag_p_l)/wl[0];
  enzo_float h_r = (Ur[4] + wr[4] + mag_p_r)/wr[0];
  enzo_float h_roe = (sqrtrho_l * h_l + sqrtrho_r * h_r) * coef;

  // Magnetic Fields (formulas are different pattern from above)
  enzo_float bi_roe = wl[5]; // equal to wr[5]
  enzo_float bj_roe = (sqrtrho_l * wr[6] + sqrtrho_r * wl[6])*coef;
  enzo_float bk_roe = (sqrtrho_l * wr[7] + sqrtrho_r * wl[7])*coef;

  
  // Calculate fast magnetosonic speed for Roe-averaged quantities (eqn B18)
  enzo_float gamma_prime = eos->get_gamma()-1.;
  enzo_float x_prime = ((std::pow(wl[6]-wr[6],2) + std::pow(wl[7]-wr[7],2) )
			* 0.5 * (gamma_prime-1) * coef);
  enzo_float y_prime = (gamma_prime-1)*(wl[0]+wr[0])*0.5/rho_roe;

  enzo_float v_roe2 = vi_roe*vi_roe + vj_roe*vj_roe + vk_roe*vk_roe;
  enzo_float b_roe2 = bi_roe*bi_roe + bj_roe*bj_roe + bk_roe*bk_roe;

  // Need tidle_a^2
  enzo_float tilde_a2 = (gamma_prime * (h_roe - 0.5* v_roe2 - b_roe2/rho_roe)
			 - x_prime);

  // Need Roe averaged square Alfven speed (along ith dim and total magnitude)
  enzo_float tilde_vai2 = bi_roe * bi_roe / rho_roe;
  enzo_float tilde_va2 = (tilde_vai2 + (gamma_prime - y_prime) *
			  (bj_roe * bj_roe + bk_roe * bk_roe));

  enzo_float cfast = std::sqrt(0.5 * (tilde_a2 + tilde_va2)
			       + std::sqrt(std::pow(tilde_a2 + tilde_va2, 2)
					   - 4 * tilde_a2 * tilde_vai2));

  *bp = std::fmax(std::fmax(vi_roe + cfast, right_speed), 0.);
  *bm = std::fmin(std::fmin(vi_roe - cfast, left_speed),  0.);
}



void EnzoRiemannHLLE::interface_flux_ (enzo_float *prim, enzo_float *cons,
				       enzo_float *fluxes,
				       enzo_float mag_pressure)
{
  // can simplify things, like in Athena++ (described in Toro)
  // they use a modified interface flux form (it involves fewer steps)

  // This assumes that MHD is included
  // This may be better handled by the EquationOfState
  enzo_float rho, vi, vj, vk, p, Bi, Bj, Bk, B2, etot, momentum_i;
  rho = prim[0];
  vi = prim[1];
  vj = prim[2];
  vk = prim[3];
  p  = prim[4];
  Bi = prim[5];
  Bj = prim[6];
  Bk = prim[7];

  B2 = Bi*Bi + Bj*Bj + Bk*Bk;

  momentum_i = rho*vi;
  // Compute Fluxes
  fluxes[0] = cons[1];

  // Fluxes for Mx, My, Mz
  fluxes[1] = cons[1]*vi - Bi*Bi + p + mag_pressure;
  fluxes[2] = cons[2]*vi - Bj*Bi;
  fluxes[3] = cons[3]*vi - Bk*Bi;

  // Flux for etot
  etot = cons[4];
  fluxes[4] = (etot + p + mag_pressure)*vi - (Bi*vi + Bj*vj + Bk*vk)*Bi;

  // Fluxes for Bi,Bj,Bk
  fluxes[5] = 0;
  fluxes[6] = Bj*vi - Bi*vj;
  fluxes[7] = Bk*vi - Bi*vk;
  
}
