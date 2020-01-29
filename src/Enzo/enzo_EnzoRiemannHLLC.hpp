// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannHLLC.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon Nov 25 2019
/// @brief    [\ref Enzo] Enzo's HLLC approximate Riemann Solver. Adapted from
/// the original Enzo's flux_hllc.F, written by John Wise. Some of the
/// formatting is modelled after EnzoRiemannHLLE. It does not include diffusion
/// terms.

/// We are intentionally skipping the diffusion terms. These require the use of
/// cell-centered values and can be split into a separate step

#ifndef ENZO_ENZO_RIEMANN_HLLC_HPP
#define ENZO_ENZO_RIEMANN_HLLC_HPP

struct HLLCImpl
{
  /// @class    EnzoRiemannHLLC
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of HLLC approximate Riemann
  /// Solver. This should not be used with isothermal equations of state.
public:
  static int scratch_space_length(const int n_cons_keys) { }

  static void calc_riemann_fluxes
  (const enzo_float flux_l[], const enzo_float flux_r[],
   const enzo_float prim_l[], const enzo_float prim_r[],
   const enzo_float cons_l[], const enzo_float cons_r[],
   const enzo_float pressure_l, const enzo_float pressure_r,
   const EnzoAdvectionFieldLUT lut, const int n_keys,
   const bool barotropic_eos, const enzo_float gamma,
   const enzo_float isothermal_cs, const bool dual_energy,
   const int iz, const int iy, const int ix, EFlt3DArray flux_arrays[],
   enzo_float scratch_space[], enzo_float &vi_bar) throw()
  {

    ASSERT("HLLCImpl::calc_riemann_fluxes",
	   "HLLC should not be used for barotropic fluids",
	   !barotropic_eos);
    
    // THE FOLLOWING SHOULD BE FACTORED OUT INTO ITS OWN FUNCTION THAT GETS
    // SHARED BY HLLE (when used WITHOUT MHD):

    enzo_float gamma1 = gamma - 1.;

    // compute min max eigenvalues
    //     - per eqn B17 these are Roe averaged velocity in the ith direction
    //       minus/plus Roe averaged fast magnetosonic wavespeed
    //     - Will probably offload this to a method of EnzoEquationOfState
    enzo_float sqrtrho_l = std::sqrt(prim_l[lut.density]);
    enzo_float sqrtrho_r = std::sqrt(prim_r[lut.density]);
    enzo_float coef = 1.0/(sqrtrho_l + sqrtrho_r);

    // density and velocity
    enzo_float rho_roe = sqrtrho_l*sqrtrho_r;
    enzo_float vi_roe = (sqrtrho_l * prim_l[lut.velocity_i] +
			 sqrtrho_r * prim_r[lut.velocity_i])*coef;
    enzo_float vj_roe = (sqrtrho_l * prim_l[lut.velocity_j] +
			 sqrtrho_r * prim_r[lut.velocity_j])*coef;
    enzo_float vk_roe = (sqrtrho_l * prim_l[lut.velocity_k] +
			 sqrtrho_r * prim_r[lut.velocity_k])*coef;

    enzo_float v_roe2 = vi_roe*vi_roe + vj_roe*vj_roe + vk_roe*vk_roe;


    // The enthalpy H=(E+P)/d is averaged for adiabatic flows rather than
    // averaging E or P directly
    // sqrtql*hl = sqrtrho_l*(el+pl)/dl = (el+pl)/sqrtrho_l

    //etot_dens_l = cons_l[lut.total_energy];
    //etot_dens_r = cons_r[lut.total_energy];
    //h_roe = ((etot_dens_l + pressure_l)/sqrtrho_l + 
    //	     (etot_dens_r + pressure_r)/sqrtrho_r) * coef

    // enthalpy:  h = (etot + thermal pressure + magnetic pressure) / rho
    // here etot is total energy density.
    // (we could refactor the following to take advantage of the fact that
    //  prim_l[lut.total_energy] already has density divided out)
    enzo_float h_l = ((cons_l[lut.total_energy] + pressure_l) /
		      prim_l[lut.density]);
    enzo_float h_r = ((cons_r[lut.total_energy] + pressure_r) /
		      prim_r[lut.density]);
    enzo_float h_roe = (sqrtrho_l * h_l + sqrtrho_r * h_r) * coef;

    // Compute characteristics using Roe-averaged values:
    enzo_float tiny2 = 0.; // tiny number to serve as a floor
    enzo_float cs = std::sqrt(gamma1 * std::max(h_roe - 0.5 * v_roe2, tiny2));
    // char1 = vi_roe - cs;       char2 = vi_roe + cs;

    // compute min/max speeds:

    enzo_float cs_l0 = EnzoRiemann::sound_speed_(prim_l, lut, pressure_l, gamma);
    enzo_float cs_r0 = EnzoRiemann::sound_speed_(prim_r, lut, pressure_r, gamma);
    enzo_float left_speed = (prim_l[lut.velocity_i] - cs_l0);
    enzo_float right_speed =  (prim_r[lut.velocity_i] + cs_r0);

    enzo_float cs_l = std::min(left_speed, vi_roe - cs);
    enzo_float cs_r = std::max(right_speed, vi_roe + cs);

    enzo_float bm = std::fmin(cs_l, 0.0);
    enzo_float bp = std::fmax(cs_r, 0.0);

    // THIS IS WHERE THE DUPLICATION ENDS

    // Compute the contact wave speed (cw) and pressure (cp).  We do this
    // for all cases because we need to correct the momentum and energy
    // flux along the contact.

    enzo_float tl = pressure_l - (cs_l - prim_l[lut.velocity_i]) * prim_l[lut.density] * prim_l[lut.velocity_i];
    enzo_float tr = pressure_r - (cs_r - prim_r[lut.velocity_i]) * prim_r[lut.density] * prim_r[lut.velocity_i];
    enzo_float dl =  prim_l[lut.density] * (cs_l - prim_l[lut.velocity_i]);
    enzo_float dr = -prim_r[lut.density] * (cs_r - prim_r[lut.velocity_i]);
    enzo_float q1 = 1.0 / (dl+dr);
    // cw is given by Toro (10.37)
    enzo_float cw = (tr - tl)*q1;
    // the following is rearranged from Toro (10.42)
    enzo_float cp = (dl*tr + dr*tl)*q1;


    // Compute the weights for fluxes
    enzo_float sl, sr, sm;
    if (cw >= 0.) {
      sl =  cw / (cw - bm);
      sr = 0.;
      sm = -bm / (cw - bm);
    } else {
      sl = 0.;
      sr = -cw / (bp - cw);
      sm =  bp / (bp - cw);
    }


    // apply floor to contact pressure
    cp = std::max(cp, 0.);
    // The following was commented out in the original code:
    // if ((sm > 0.0) && (cp < 0.0)){
    //   WARNING4("HLLCImpl::calc_riemann_fluxes",
    //            "Negative contact pressure at (iz=%d,iy=%d,ix=%d): %e",
    //            iz, iy, ix, cp);
    //   cp = 0.0;
    // }


    // Compute the left and right fluxes along the characteristics.
    // For dual energy formalism, treat internal energy like an advected
    // quantity (source term is handled separately)

    enzo_float momentumi_l = cons_l[lut.velocity_i];
    enzo_float momentumi_r = cons_r[lut.velocity_i];

    enzo_float dfl,dfr, ufl,ufr, vfl,vfr, wfl,wfr, efl,efr;

    dfl = momentumi_l - bm*prim_l[lut.density];
    dfr = momentumi_r - bp*prim_r[lut.density];

    ufl = momentumi_l * (prim_l[lut.velocity_i] - bm) + pressure_l;
    ufr = momentumi_r * (prim_r[lut.velocity_i] - bp) + pressure_r;

    vfl = prim_l[lut.density]*prim_l[lut.velocity_j] * (prim_l[lut.velocity_i] - bm);
    vfr = prim_r[lut.density]*prim_r[lut.velocity_j] * (prim_r[lut.velocity_i] - bp);

    wfl = prim_l[lut.density]*prim_l[lut.velocity_k] * (prim_l[lut.velocity_i] - bm);
    wfr = prim_r[lut.density]*prim_r[lut.velocity_k] * (prim_r[lut.velocity_i] - bp);

    efl = (cons_l[lut.total_energy] * (prim_l[lut.velocity_i] - bm)
	   + pressure_l*prim_l[lut.velocity_i]);
    efr = (cons_r[lut.total_energy] * (prim_r[lut.velocity_i] - bp)
	   + pressure_r*prim_r[lut.velocity_i]);

    // internal energy (Gas energy) is treated like a passively advected
    // quantity (like transverse velocity components and color fields)
    enzo_float eint_fl, eint_fr;
    if (dual_energy){
       eint_fl = ((prim_l[lut.velocity_i]-bm) * prim_l[lut.internal_energy]
		  * prim_l[lut.density]);
       eint_fr = ((prim_r[lut.velocity_i]-bp) * prim_r[lut.internal_energy]
		  * prim_r[lut.density]);
    }
    // An aside: compute the interface velocity (that might be used to compute
    // the internal energy source term)

    vi_bar = (sl * (prim_l[lut.velocity_i] - bm) +
	      sr * (prim_r[lut.velocity_i] - bp));

    // compute HLLC Flux at interface (without diffusion)
    flux_arrays[lut.density](iz,iy,ix) = sl*dfl + sr*dfr;
    flux_arrays[lut.velocity_i](iz,iy,ix) = sl*ufl + sr*ufr;
    flux_arrays[lut.velocity_j](iz,iy,ix) = sl*vfl + sr*vfr;
    flux_arrays[lut.velocity_k](iz,iy,ix) = sl*wfl + sr*wfr;
    flux_arrays[lut.total_energy](iz,iy,ix) = sl*efl + sr*efr;

    if (dual_energy){
      flux_arrays[lut.internal_energy](iz,iy,ix) = sl*eint_fl + sr*eint_fr;
    }

    // Add the weighted contribution of the flux along the contact for
    // velocity_i and total_energy
    flux_arrays[lut.velocity_i](iz,iy,ix) += (sm * cp);
    flux_arrays[lut.total_energy](iz,iy,ix) += (sm * cp * cw);

    flux_arrays[lut.bfield_i](iz,iy,ix) = 0;
    flux_arrays[lut.bfield_j](iz,iy,ix) = 0;
    flux_arrays[lut.bfield_k](iz,iy,ix) = 0;
  }

};

/// @class    EnzoRiemannHLLC
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLC approximate Riemann Solver
using EnzoRiemannHLLC = EnzoRiemannImpl<HLLCImpl>;

#endif /* ENZO_ENZO_RIEMANN_HLLC_HPP */
