// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannHLLE.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of HLLE approximate Riemann Solver

#ifndef ENZO_ENZO_RIEMANN_HLLE_HPP
#define ENZO_ENZO_RIEMANN_HLLE_HPP

#include <cfloat>
#include <cmath>


// Implementation Notes:
//   - Based off of the HLLE Riemann Solver described in Stone+08
//   - Currently only supports adiabatic fluids
//   - This class is primarily implemented to allow for easy extension to
//     incorporate CRs (just subclass and overwrite the wave_speeds_ method)
// Misc Notes:
//   - Stone+ (08) has a typo in eqn 52: (q_i-q_{i-1}) should be (q_R-q_L)
//     The modified formula is actually used in Athena++


struct AthenaHLLEWavespeed
{
  /// @class    AthenaHLLEWavespeed
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the wavespeed calculation described by
  ///           Stone+08

  void operator()(const enzo_float wl[], const enzo_float wr[],
		  const enzo_float Ul[], const enzo_float Ur[],
		  const field_lut prim_lut, const field_lut cons_lut,
		  const enzo_float gamma, enzo_float *bp, enzo_float *bm)
  {
    field_lut plut = prim_lut;
    field_lut clut = cons_lut;
    // Calculate wavespeeds as specified by S4.3.1 of Stone+08
    // if LM and L0 are max/min eigenvalues of Roe's matrix:
    //       bp = max(LM, vi_r + cf_r, 0) (eqn 53)
    //       bm = max(L0, vi_l - cf_l, 0) (eqn 54)
    // vi_r and vi_l are velocities of left and right states along ith dimension
    // cf_r and cf_l are the fast magnetosonic speeds in the left and right
    // states

    // First, compute left and right speeds
    //    left_speed = vi_l - cf_l;      right_speed = vi_r + cf_r
    enzo_float left_speed = (wl[plut.velocity_i] -
			     EnzoRiemann::fast_magnetosonic_speed_(wl, plut,
								   gamma));
    enzo_float right_speed = (wr[plut.velocity_i] +
			      EnzoRiemann::fast_magnetosonic_speed_(wr, plut,
								    gamma));

    // Next, compute min max eigenvalues
    //     - per eqn B17 these are Roe averaged velocity in the ith direction
    //       minus/plus Roe averaged fast magnetosonic wavespeed
    //     - Will probably offload this to a method of EnzoEquationOfState
    enzo_float sqrtrho_l = std::sqrt(wl[plut.density]);
    enzo_float sqrtrho_r = std::sqrt(wr[plut.density]);
    enzo_float coef = 1.0/(sqrtrho_l + sqrtrho_r);

    // density and velocity
    enzo_float rho_roe = sqrtrho_l*sqrtrho_r;
    enzo_float vi_roe = (sqrtrho_l * wl[plut.velocity_i] +
			 sqrtrho_r * wr[plut.velocity_i])*coef;
    enzo_float vj_roe = (sqrtrho_l * wl[plut.velocity_j] +
			 sqrtrho_r * wr[plut.velocity_j])*coef;
    enzo_float vk_roe = (sqrtrho_l * wl[plut.velocity_k] +
			 sqrtrho_r * wr[plut.velocity_k])*coef;

    enzo_float mag_p_l = EnzoRiemann::mag_pressure_(wl, plut);
    enzo_float mag_p_r = EnzoRiemann::mag_pressure_(wr, plut);

    // enthalpy:  h = (etot + thermal pressure + magnetic pressure) / rho
    enzo_float h_l = ((Ul[clut.total_energy] + wl[plut.pressure] + mag_p_l) /
		      wl[plut.density]);
    enzo_float h_r = ((Ur[clut.total_energy] + wr[plut.pressure] + mag_p_r) /
		      wr[plut.density]);
    enzo_float h_roe = (sqrtrho_l * h_l + sqrtrho_r * h_r) * coef;

    // Magnetic Fields (formulas are different pattern from above)
    enzo_float bi_roe = wl[plut.bfield_i]; // equal to wr["bfield_i"]
    enzo_float bj_roe = (sqrtrho_l * wr[plut.bfield_j] +
			 sqrtrho_r * wl[plut.bfield_j])*coef;
    enzo_float bk_roe = (sqrtrho_l * wr[plut.bfield_k] +
			 sqrtrho_r * wl[plut.bfield_k])*coef;


    // Calculate fast magnetosonic speed for Roe-averaged quantities (eqn B18)
    enzo_float gamma_prime = gamma - 1.;
    enzo_float x_prime = ((std::pow(wl[plut.bfield_j]-wr[plut.bfield_j],2) +
			   std::pow(wl[plut.bfield_k]-wr[plut.bfield_k],2) )
			  * 0.5 * (gamma_prime-1) * coef);
    enzo_float y_prime = ((gamma_prime-1)*(wl[plut.density]+wr[plut.density])
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
    //CkPrintf("left_speed = %.15lf, right_speed = %.15lf\n",left_speed, right_speed);
    //fflush(stdout);

    *bp = std::fmax(std::fmax(vi_roe + cfast, right_speed), 0.);
    *bm = std::fmin(std::fmin(vi_roe - cfast, left_speed),  0.);
  }
};


template <class WaveSpeedFunctor>
struct HLLEImpl
{
  /// @class    EnzoRiemannHLLE
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of HLLE approximate Riemann
  /// Solver. 
public:
  static int scratch_space_length(const int n_cons_keys){
    return 0;
  }

  static void calc_riemann_fluxes
  (const enzo_float flux_l[], const enzo_float flux_r[],
   const enzo_float prim_l[], const enzo_float prim_r[],
   const enzo_float cons_l[], const enzo_float cons_r[],
   const field_lut prim_lut, const field_lut cons_lut, const int n_keys,
   const bool barotropic_eos, const enzo_float gamma,
   const enzo_float isothermal_cs, const int iz, const int iy, const int ix,
   EFlt3DArray flux_arrays[], enzo_float scratch_space[])
  { 
    // there is no scratch_space
    WaveSpeedFunctor wave_speeds;
    //AthenaHLLEWavespeed wave_speeds;
    enzo_float bp, bm;

    // Compute wave speeds
    wave_speeds(prim_l, prim_r, cons_l, cons_r, prim_lut, cons_lut, gamma,
		&bp, &bm);

    // Compute the actual riemann fluxes
    for (int field = 0; field<n_keys; field++){
      flux_arrays[field](iz,iy,ix) = ((bp*flux_l[field] - bm*flux_r[field] +
				       (cons_r[field] - cons_l[field])*bp*bm)
				      / (bp - bm));
    }
  }
};


/// @class    EnzoRiemannHLLE
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLE approximate Riemann Solver
using EnzoRiemannHLLE = EnzoRiemannImpl<HLLEImpl<AthenaHLLEWavespeed>>;

// To define a Riemann Solver for cosmic ray transport, we would probably
// define a separate WaveSpeedFunctor struct and define a separate HLLE
// Riemann Solver (maybe EnzoRiemannCRHLLE)

#endif /* ENZO_ENZO_RIEMANN_HLLE_HPP */
