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

struct HLLCKernel
{
  /// @class    HLLCKernel
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of HLLC approximate Riemann
  /// Solver. This should not be used with isothermal equations of state.

public: // typedefs
  using WaveSpeedFunctor = EinfeldtWavespeed<HydroLUT>;
  using LUT = typename WaveSpeedFunctor::LUT;
  using EOSStructT = EnzoEOSIdeal;

public: // fields
  const KernelConfig<EOSStructT> config;

public: // methods

  FORCE_INLINE void operator()(const int iz,
                               const int iy,
                               const int ix) const noexcept
  {

    const int external_velocity_i = config.dim + LUT::velocity_i;
    const int external_velocity_j = ((config.dim+1)%3) + LUT::velocity_i;
    const int external_velocity_k = ((config.dim+2)%3) + LUT::velocity_i;

    // load primitives into local variables
    lutarray<LUT> prim_l;
    prim_l[LUT::density]     = config.prim_arr_l(LUT::density,iz,iy,ix);
    prim_l[LUT::velocity_i]  = config.prim_arr_l(external_velocity_i,iz,iy,ix);
    prim_l[LUT::velocity_j]  = config.prim_arr_l(external_velocity_j,iz,iy,ix);
    prim_l[LUT::velocity_k]  = config.prim_arr_l(external_velocity_k,iz,iy,ix);
    // this actually stores pressure:
    prim_l[LUT::total_energy]= config.prim_arr_l(LUT::total_energy,iz,iy,ix);

    lutarray<LUT> prim_r;
    prim_r[LUT::density]     = config.prim_arr_r(LUT::density,iz,iy,ix);
    prim_r[LUT::velocity_i]  = config.prim_arr_r(external_velocity_i,iz,iy,ix);
    prim_r[LUT::velocity_j]  = config.prim_arr_r(external_velocity_j,iz,iy,ix);
    prim_r[LUT::velocity_k]  = config.prim_arr_r(external_velocity_k,iz,iy,ix);
    // this actually stores pressure:
    prim_r[LUT::total_energy]= config.prim_arr_r(LUT::total_energy,iz,iy,ix);

    // load left and right pressure values
    const enzo_float pressure_l = prim_l[LUT::total_energy];
    const enzo_float pressure_r = prim_r[LUT::total_energy];

    // get the conserved quantities
    const lutarray<LUT> cons_l
      = enzo_riemann_utils::compute_conserved<LUT>(prim_l, config.eos);
    const lutarray<LUT> cons_r
      = enzo_riemann_utils::compute_conserved<LUT>(prim_r, config.eos);

    enzo_float cs_l,cs_r;
    WaveSpeedFunctor wave_speeds;
    wave_speeds(prim_l, prim_r, cons_l, cons_r, pressure_l, pressure_r,
		config.eos, &cs_r, &cs_l);

    enzo_float bm = std::fmin(cs_l, 0.0);
    enzo_float bp = std::fmax(cs_r, 0.0);

    // Compute the contact wave speed (cw) and pressure (cp).  We do this
    // for all cases because we need to correct the momentum and energy
    // flux along the contact.

    enzo_float tl = (pressure_l - (cs_l - prim_l[LUT::velocity_i]) *
                     prim_l[LUT::density] * prim_l[LUT::velocity_i]);
    enzo_float tr = (pressure_r - (cs_r - prim_r[LUT::velocity_i]) *
                     prim_r[LUT::density] * prim_r[LUT::velocity_i]);
    enzo_float dl =  prim_l[LUT::density] * (cs_l - prim_l[LUT::velocity_i]);
    enzo_float dr = -prim_r[LUT::density] * (cs_r - prim_r[LUT::velocity_i]);
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
    cp = std::max(cp, (enzo_float)0.);
    // The following was commented out in the original code:
    // if ((sm > 0.0) && (cp < 0.0)){
    //   WARNING4("HLLCImpl::calc_riemann_fluxes",
    //            "Negative contact pressure at (iz=%d,iy=%d,ix=%d): %e",
    //            iz, iy, ix, cp);
    //   cp = 0.0;
    // }


    // Compute the left and right fluxes along the characteristics.

    enzo_float momentumi_l = cons_l[LUT::velocity_i];
    enzo_float momentumi_r = cons_r[LUT::velocity_i];

    enzo_float dfl,dfr, ufl,ufr, vfl,vfr, wfl,wfr, efl,efr;

    dfl = momentumi_l - bm*prim_l[LUT::density];
    dfr = momentumi_r - bp*prim_r[LUT::density];

    ufl = momentumi_l * (prim_l[LUT::velocity_i] - bm) + pressure_l;
    ufr = momentumi_r * (prim_r[LUT::velocity_i] - bp) + pressure_r;

    vfl = (prim_l[LUT::density] * prim_l[LUT::velocity_j] *
           (prim_l[LUT::velocity_i] - bm));
    vfr = (prim_r[LUT::density] * prim_r[LUT::velocity_j] *
           (prim_r[LUT::velocity_i] - bp));

    wfl = (prim_l[LUT::density] * prim_l[LUT::velocity_k] *
           (prim_l[LUT::velocity_i] - bm));
    wfr = (prim_r[LUT::density] * prim_r[LUT::velocity_k] *
           (prim_r[LUT::velocity_i] - bp));

    efl = (cons_l[LUT::total_energy] * (prim_l[LUT::velocity_i] - bm)
	   + pressure_l * prim_l[LUT::velocity_i]);
    efr = (cons_r[LUT::total_energy] * (prim_r[LUT::velocity_i] - bp)
	   + pressure_r * prim_r[LUT::velocity_i]);

    // compute HLLC Flux at interface (without diffusion)
    config.flux_arr(LUT::density,iz,iy,ix) = sl*dfl + sr*dfr;
    config.flux_arr(external_velocity_i,iz,iy,ix) = sl*ufl + sr*ufr;
    config.flux_arr(external_velocity_j,iz,iy,ix) = sl*vfl + sr*vfr;
    config.flux_arr(external_velocity_k,iz,iy,ix) = sl*wfl + sr*wfr;
    config.flux_arr(LUT::total_energy,iz,iy,ix) = sl*efl + sr*efr;

    // Add the weighted contribution of the flux along the contact for
    // velocity_i and total_energy
    // (if you break 10.44 into 2 fractions, we are adding the right one)
    config.flux_arr(external_velocity_i,iz,iy,ix) += (sm * cp);
    config.flux_arr(LUT::total_energy,iz,iy,ix) += (sm * cp * cw);


    // finally, deal with dual energy stuff:
    // compute passive advection internal energy flux
    config.internal_energy_flux_arr(iz,iy,ix) =
      enzo_riemann_utils::passive_eint_flux
      (prim_l[LUT::density], pressure_l, prim_r[LUT::density], pressure_r,
       config.eos, config.flux_arr(LUT::density,iz,iy,ix));

    // An aside: compute the interface velocity (this is used to compute the
    // internal energy source term)
    config.velocity_i_bar_arr(iz,iy,ix) =
      (sl * (prim_l[LUT::velocity_i] - bm) +
       sr * (prim_r[LUT::velocity_i] - bp));

  }

};

/// @class    EnzoRiemannHLLC
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLC approximate Riemann Solver
using EnzoRiemannHLLC = EnzoRiemannImpl<HLLCKernel>;

#endif /* ENZO_ENZO_RIEMANN_HLLC_HPP */
