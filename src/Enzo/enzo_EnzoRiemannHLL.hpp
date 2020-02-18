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

template<bool MHD>
struct EinfeldtWavespeed
{
  /// @class    EinfeldtWavespeed
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the wavespeed calculation described by
  ///           Stone+08. Toro attributes this scheme to Einfeldt 1988. When
  ///           used for the HLL solver, it makes the HLLE solver. This scheme
  ///           allows the min/max eigenvalues of Roe's matrix to be
  ///           wavespeeds.
  ///
  /// This estimator provides wave speeds that bound almost all possible
  /// physical signal speeds. According to Batten, Clarke, Lambert, & Causon
  /// (1997), this doesn't account for the case when either of the outermost
  /// waves in the Riemann Fan is a shock. Nevertheless, they find this is a
  /// robust method and recommend it in most cases

  void operator()(const enzo_float wl[], const enzo_float wr[],
		  const enzo_float Ul[], const enzo_float Ur[],
		  const enzo_float pressure_l, const enzo_float pressure_r,
		  const EnzoAdvectionFieldLUT lut,  const enzo_float gamma,
		  enzo_float *bp, enzo_float *bm) throw()
  {
    // Calculate wavespeeds as specified by S4.3.1 of Stone+08
    // if LM and L0 are max/min eigenvalues of Roe's matrix:
    //       bp = max(LM, vi_r + c_r, 0) (eqn 53)
    //       bm = max(L0, vi_l - c_l, 0) (eqn 54)
    // vi_r and vi_l are velocities of left and right states along ith dimension
    // c_r and c_l are the maximum wavespeeds (sound speeds if hydro, fast
    // magnetosonic speed if MHD) in the left and right states

    // First, compute left and right speeds
    //    left_speed = vi_l - c_l;      right_speed = vi_r + c_r
    enzo_float c_l, c_r;
    if (MHD){
      c_l = EnzoRiemann::fast_magnetosonic_speed_(wl, lut, pressure_l, gamma);
      c_r = EnzoRiemann::fast_magnetosonic_speed_(wr, lut, pressure_r, gamma);
    } else {
      c_l = EnzoRiemann::sound_speed_(wl, lut, pressure_l, gamma);
      c_r = EnzoRiemann::sound_speed_(wr, lut, pressure_r, gamma);
    }
    enzo_float left_speed = (wl[lut.velocity_i] - c_l);
    enzo_float right_speed = (wr[lut.velocity_i] + c_r);

    // Next, compute min, max eigenvalues of Roe matrix. Per eqns B2 and B17
    // these are the Roe averaged velocity in the ith direction minus/plus
    // the maximum wavespeeds Roe averaged signal speed (sound speed for hydro
    // and fast magnetosonic speed for MHD)
    enzo_float sqrtrho_l = std::sqrt(wl[lut.density]);
    enzo_float sqrtrho_r = std::sqrt(wr[lut.density]);
    enzo_float inv_sqrtrho_tot = 1.0/(sqrtrho_l + sqrtrho_r);

    // Compute Roe-averaged velocity:
    enzo_float vi_roe = (sqrtrho_l * wl[lut.velocity_i] +
			 sqrtrho_r * wr[lut.velocity_i]) * inv_sqrtrho_tot;
    enzo_float vj_roe = (sqrtrho_l * wl[lut.velocity_j] +
			 sqrtrho_r * wr[lut.velocity_j]) * inv_sqrtrho_tot;
    enzo_float vk_roe = (sqrtrho_l * wl[lut.velocity_k] +
			 sqrtrho_r * wr[lut.velocity_k]) * inv_sqrtrho_tot;
    enzo_float v_roe2 = vi_roe*vi_roe + vj_roe*vj_roe + vk_roe*vk_roe;

    // Compute Roe-averaged specific enthalpy
    enzo_float ptot_l = pressure_l;
    enzo_float ptot_r = pressure_r;

    if (MHD){
      ptot_l += EnzoRiemann::mag_pressure_(wl, lut);
      ptot_r += EnzoRiemann::mag_pressure_(wr, lut);
    }

    // enthalpy:  h = (etot + thermal pressure + magnetic pressure) / rho
    // here etot is total energy density.
    // (we could refactor the following to take advantage of the fact that
    //  wl[lut.total_energy] already has density divided out)
    enzo_float h_l = (Ul[lut.total_energy] + ptot_l) / wl[lut.density];
    enzo_float h_r = (Ur[lut.total_energy] + ptot_r) / wr[lut.density];
    enzo_float h_roe = (sqrtrho_l * h_l + sqrtrho_r * h_r) * inv_sqrtrho_tot;

    enzo_float c_roe;
    if (MHD){
      c_roe = roe_cfast_(wl, wr, sqrtrho_l, sqrtrho_r, inv_sqrtrho_tot, v_roe2,
			 h_roe, gamma, lut);
    } else {
      c_roe = roe_cs_(v_roe2, h_roe, gamma);
    }

    *bp = std::fmax(vi_roe + c_roe, right_speed);
    *bm = std::fmin(vi_roe - c_roe, left_speed);
  }

  /// Return the sound speed needed to get the max and min eigenvalues for
  /// Roe's matrix for adiabatic Hydrodynamics
  ///
  /// @param v_roe2 The squared magnitude of the roe averaged velocity
  /// @param h_roe The Roe averaged specific enthalpy.
  /// @param gamma The adiabatic index
  enzo_float roe_cs_(const enzo_float v_roe2, const enzo_float h_roe,
		     const enzo_float gamma) const throw()
  {
    enzo_float tiny2 = 0.;
    return std::sqrt((gamma - 1) * std::max(h_roe - 0.5 * v_roe2, tiny2));
  }

  /// Return the fast magnetosonic speed needed to get the max and min
  /// eigenvalues for Roe's matrix for adiabatic Magnetohydrodynamics
  enzo_float roe_cfast_(const enzo_float wl[], const enzo_float wr[],
			const enzo_float sqrtrho_l, const enzo_float sqrtrho_r,
			const enzo_float inv_sqrtrho_tot,
			const enzo_float v_roe2, const enzo_float h_roe,
			const enzo_float gamma,
			const EnzoAdvectionFieldLUT lut) const throw()
  {
    // Roe-averaged density & B-field (formulas different velocity & enthalpy)
    enzo_float rho_roe = sqrtrho_l*sqrtrho_r;
    enzo_float bi_roe = wl[lut.bfield_i]; // equal to wr["bfield_i"]
    enzo_float bj_roe = (sqrtrho_l * wr[lut.bfield_j] +
			 sqrtrho_r * wl[lut.bfield_j]) * inv_sqrtrho_tot;
    enzo_float bk_roe = (sqrtrho_l * wr[lut.bfield_k] +
			 sqrtrho_r * wl[lut.bfield_k]) * inv_sqrtrho_tot;
    enzo_float b_roe2 = bi_roe*bi_roe + bj_roe*bj_roe + bk_roe*bk_roe;

    // Calculate fast magnetosonic speed for Roe-averaged quantities (eqn B18)
    enzo_float gamma_prime = gamma - 1.;
    enzo_float x_prime = ((std::pow(wl[lut.bfield_j]-wr[lut.bfield_j],2) +
			   std::pow(wl[lut.bfield_k]-wr[lut.bfield_k],2) )
			  * 0.5 * (gamma_prime-1) * inv_sqrtrho_tot);
    enzo_float y_prime = ((gamma_prime-1)*(wl[lut.density]+wr[lut.density])
			  *0.5/rho_roe);

    // Need tidle_a^2
    enzo_float tilde_a2 = (gamma_prime * (h_roe - 0.5* v_roe2 - b_roe2/rho_roe)
			   - x_prime);

    // Need Roe averaged square Alfven speed (along ith dim and total magnitude)
    enzo_float tilde_vai2 = bi_roe * bi_roe / rho_roe;
    enzo_float tilde_va2 = (tilde_vai2 + (gamma_prime - y_prime) *
			    (bj_roe * bj_roe + bk_roe * bk_roe) / rho_roe);

    return std::sqrt(0.5 * (tilde_a2 + tilde_va2 +
			    std::sqrt(std::pow(tilde_a2 + tilde_va2, 2) -
				      4 * tilde_a2 * tilde_vai2)));
  }
};


struct DavisWavespeedMHD
{
  /// @class    DavisHLLEWavespeedMHD
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates one of the wavespeed estimators
  ///           described by Davis 1998. It was used in Enzo's Runge-Kutta MHD
  ///           HLL Riemann solver in Riemann_HLL_MHD.C which was (written by
  ///           Peng Wang)
  ///
  /// Toro, describes this version as the "most well-known approach." He goes
  /// on to describe it as "exceedling simple" and not recommended for
  /// practical calculations. He attributes this to the scheme from Davis 1988
public:
  void operator()(const enzo_float wl[], const enzo_float wr[],
		  const enzo_float Ul[], const enzo_float Ur[],
		  const enzo_float pressure_l, const enzo_float pressure_r,
		  const EnzoAdvectionFieldLUT lut,  const enzo_float gamma,
		  enzo_float *bp, enzo_float *bm) throw()
  {
    // Should be relatively easy to adapt and make compatible with barotropic
    // eos
    ERROR("EnzoHLLEWavespeed","This hasn't been tested yet");

    // compute left and right fast magnetosonic speed assuming that cos2 = 1
    // not sure why cos2 = 1 was chosen, (cos2 = 0, would yield faster speeds)
    enzo_float cf_l, cf_r;
    cf_l = EnzoRiemann::fast_magnetosonic_speed_(wl, lut, pressure_l, gamma, 1);
    cf_r = EnzoRiemann::fast_magnetosonic_speed_(wr, lut, pressure_r, gamma, 1);

    enzo_float lp_l = wl[lut.velocity_i] + cf_l;
    enzo_float lm_l = wl[lut.velocity_i] - cf_l;
    enzo_float lp_r = wr[lut.velocity_i] + cf_r;
    enzo_float lm_r = wr[lut.velocity_i] - cf_r;

    // The following is equivalent to the enzo version except bm has been
    // multiplied by -1
    *bp = std::fmax(lp_l, lp_r);
    *bm = std::fmin(lm_l, lm_r);
  }
};


template <class WaveSpeedFunctor>
struct HLLImpl
{
  /// @class    EnzoRiemannHLLE
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of HLLE approximate Riemann
  /// Solver. 
public:
  static int scratch_space_length(const int n_cons_keys) throw(){
    return 0;
  }

  static void calc_riemann_fluxes
  (const enzo_float flux_l[], const enzo_float flux_r[],
   const enzo_float prim_l[], const enzo_float prim_r[],
   const enzo_float cons_l[], const enzo_float cons_r[],
   const enzo_float pressure_l, const enzo_float pressure_r,
   const EnzoAdvectionFieldLUT lut, const int n_keys,
   const bool barotropic_eos, const enzo_float gamma,
   const enzo_float isothermal_cs, const bool dual_energy_formalism,
   const int iz, const int iy, const int ix, EFlt3DArray flux_arrays[],
   enzo_float scratch_space[], enzo_float &vi_bar) throw()
  { 
    // there is no scratch_space
    WaveSpeedFunctor wave_speeds;
    enzo_float bp, bm;

    // Compute wave speeds
    wave_speeds(prim_l, prim_r, cons_l, cons_r, pressure_l, pressure_r, lut,
		gamma, &bp, &bm);
    bp = std::fmax(bp,0.0);
    bm = std::fmin(bm,0.0);
    enzo_float inv_speed_diff = 1./(bp - bm);

    // Compute the actual riemann fluxes (value of dual_energy_formalism is
    // irrelevant)
    for (int field = 0; field<n_keys; field++){
      flux_arrays[field](iz,iy,ix) = ((bp*flux_l[field] - bm*flux_r[field] +
				       (cons_r[field] - cons_l[field])*bp*bm)
				      * inv_speed_diff);
    }

    // Estimate the value of the vi (ith component of velocity) which is used
    // to compute the internal energy source term. The following is adopted
    // from Enzo's flux_hll.F. The logic reduces to:
    //   - when bm >=0 (left speed is positive) use vi_L or wl[lut.velocity_i]
    //   - when bp <=0 (right speed is negative) use vi_R or wr[lut.velocity_i]
    //   - otherwise linearly interpolate between vi_L and vi_R. Let the cell
    //     interface be at x=0. At some time t, the velocity is vi_L at x=t*bm
    //     and vi_R at x=t*bp. The factors of t cancel.
    vi_bar = (bp*prim_l[lut.velocity_i]
	      - bm*prim_r[lut.velocity_i]) * inv_speed_diff;
    // Alternatively, if vi is assumed to be constant in the intermediate zone
    // (which is an underlying assumption of HLLC & HLLD), it can be estimated
    // using equations for the values of rho and rho*vi in the intermediate
    // zone. This yields the equation for the contact wave speed in HLLC & HLLD
  }
};


/// @class    EnzoRiemannHLLEMHD
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLE approximate Riemann Solver (the HLL
///           solver using Einfeldt's wavespeed estimator). The same wavespeed
///           estimator is used in Athena's MHD HLLE Riemann Solver
using EnzoRiemannHLLEMHD = EnzoRiemannImpl<HLLImpl<EinfeldtWavespeed<true>>>;
using EnzoRiemannHLLE = EnzoRiemannImpl<HLLImpl<EinfeldtWavespeed<false>>>;

/// @class    EnzoRiemannHLLMHD
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLL approximate Riemann Solver using
///           the Davis MHD wavespeed estimator. The same wavespeed estimator
///           is used in Enzo's Runge-Kutta MHD HLL Riemann Solver
using EnzoRiemannHLLMHD = EnzoRiemannImpl<HLLImpl<DavisWavespeedMHD>>;

// To define a Riemann Solver for cosmic ray transport, we would probably
// define a separate WaveSpeedFunctor struct and define a separate HLL
// Riemann Solver (maybe EnzoRiemannCRHLL)

#endif /* ENZO_ENZO_RIEMANN_HLLE_HPP */
