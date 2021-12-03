// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannHLL2.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of HLLE approximate Riemann Solver

#ifndef ENZO_ENZO_RIEMANN_HLL2_HPP
#define ENZO_ENZO_RIEMANN_HLL2_HPP

template <class WaveSpeedFunctor>
struct HLLKernel
{
  /// @class    HLLImpl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of HLL approximate Riemann
  /// Solver.
public:
  const KernelConfig config;

public:

  using LUT = typename WaveSpeedFunctor::LUT;

  FORCE_INLINE void operator()(const int iz,
                               const int iy,
                               const int ix) const noexcept
  {
    const int external_velocity_i = config.dim + LUT::velocity_i;
    const int external_velocity_j = ((config.dim+1)%3) + LUT::velocity_i;
    const int external_velocity_k = ((config.dim+2)%3) + LUT::velocity_i;
    const int external_bfield_i = config.dim + LUT::bfield_i;
    const int external_bfield_j = ((config.dim+1)%3) + LUT::bfield_i;
    const int external_bfield_k = ((config.dim+2)%3) + LUT::bfield_i;

    const enzo_float gamma = config.gamma;

    // load primitives:
    auto load_prim = [=](const CelloArray<const enzo_float,4>& prim_arr)
    {
      lutarray<LUT> out;
      out[LUT::density]      = prim_arr(LUT::density,iz,iy,ix);
      out[LUT::velocity_i]   = prim_arr(external_velocity_i,iz,iy,ix);
      out[LUT::velocity_j]   = prim_arr(external_velocity_j,iz,iy,ix);
      out[LUT::velocity_k]   = prim_arr(external_velocity_k,iz,iy,ix);
      if (LUT::has_bfields){
        out[LUT::bfield_i] = prim_arr(external_bfield_i,iz,iy,ix);
        out[LUT::bfield_j] = prim_arr(external_bfield_j,iz,iy,ix);
        out[LUT::bfield_k] = prim_arr(external_bfield_k,iz,iy,ix);
      }
      // this actually stores pressure:
      out[LUT::total_energy] = prim_arr(LUT::total_energy,iz,iy,ix);
      return out;
    };

    const lutarray<LUT> prim_l = load_prim(config.prim_arr_l);
    const lutarray<LUT> prim_r = load_prim(config.prim_arr_r);

    // load left and right pressure values
    const enzo_float pressure_l = prim_l[LUT::total_energy];
    const enzo_float pressure_r = prim_r[LUT::total_energy];

    // get the conserved quantities
    const lutarray<LUT> cons_l
      = enzo_riemann_utils::compute_conserved<LUT>(prim_l, gamma);
    const lutarray<LUT> cons_r
      = enzo_riemann_utils::compute_conserved<LUT>(prim_r, gamma);

    // compute the interface fluxes
    const lutarray<LUT> flux_l
      = enzo_riemann_utils::active_fluxes<LUT>(prim_l, cons_l, pressure_l);
    const lutarray<LUT> flux_r
      = enzo_riemann_utils::active_fluxes<LUT>(prim_r, cons_r, pressure_r);

    // there is no scratch_space
    WaveSpeedFunctor wave_speeds;
    enzo_float bp, bm;

    // Compute wave speeds
    wave_speeds(prim_l, prim_r, cons_l, cons_r, pressure_l, pressure_r, gamma,
                &bp, &bm);
    bp = std::fmax(bp,0.0);
    bm = std::fmin(bm,0.0);
    enzo_float inv_speed_diff = 1./(bp - bm);

    // Compute the actual riemann fluxes
    // (we'll address dual energy considerations, afterwards)
    lutarray<LUT> flux;
    for (std::size_t field = 0; field < LUT::num_entries; field++){
      flux[field] = ((bp*flux_l[field] - bm*flux_r[field] +
                      (cons_r[field] - cons_l[field])*bp*bm) * inv_speed_diff);
    }

    config.flux_arr(LUT::density,iz,iy,ix) = flux[LUT::density];
    config.flux_arr(external_velocity_i,iz,iy,ix) = flux[LUT::velocity_i];
    config.flux_arr(external_velocity_j,iz,iy,ix) = flux[LUT::velocity_j];
    config.flux_arr(external_velocity_k,iz,iy,ix) = flux[LUT::velocity_k];
    if (LUT::has_bfields){
      config.flux_arr(external_bfield_i,iz,iy,ix) = flux[LUT::bfield_i];
      config.flux_arr(external_bfield_j,iz,iy,ix) = flux[LUT::bfield_j];
      config.flux_arr(external_bfield_k,iz,iy,ix) = flux[LUT::bfield_k];
    }
    // this actually stores pressure:
    config.flux_arr(LUT::total_energy,iz,iy,ix) = flux[LUT::total_energy];

    // finally, deal with dual energy stuff.
    // compute internal energy flux, assuming passive advection
    // (this was not handled with the rest of the fluxes)
    config.internal_energy_flux_arr(iz,iy,ix) =
      enzo_riemann_utils::passive_eint_flux
      (prim_l[LUT::density], pressure_l, prim_r[LUT::density], pressure_r,
       gamma, config.flux_arr(LUT::density,iz,iy,ix));

    // Estimate the value of the vi (ith component of velocity) which is used
    // to compute the internal energy source term. The following is adopted
    // from Enzo's flux_hll.F. The logic reduces to:
    //   - when bm>=0 (left speed is positive) use vi_L or wl[LUT::velocity_i]
    //   - when bp<=0 (right speed is negative) use vi_R or wr[LUT::velocity_i]
    //   - otherwise linearly interpolate between vi_L and vi_R. Let the cell
    //     interface be at x=0. At some time t, the velocity is vi_L at x=t*bm
    //     and vi_R at x=t*bp. The factors of t cancel.
    config.velocity_i_bar_arr(iz,iy,ix) =
      (bp*prim_l[LUT::velocity_i]
       - bm*prim_r[LUT::velocity_i]) * inv_speed_diff;
    // Alternatively, if vi is assumed to be constant in the intermediate zone
    // (which is an underlying assumption of HLLC & HLLD), it can be estimated
    // using equations for the values of rho and rho*vi in the intermediate
    // zone. This yields the equation for the contact wave speed in HLLC & HLLD
  }
};

/// @class    EnzoRiemannHLLEMHD
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLE approximate Riemann Solver (the HLL
///           solver using Einfeldt's wavespeed estimator) for MHD. The same
///           wavespeed estimator is used in Athena's MHD HLLE Riemann Solver
using EnzoRiemannHLLEMHD2 =
  EnzoRiemannImplNew<HLLKernel<EinfeldtWavespeed<MHDLUT>>>;

// @class    EnzoRiemannHLLE
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLE approximate Riemann Solver (the HLL
///           solver using Einfeldt's wavespeed estimator). The same wavespeed
///           estimator is used in Athena's MHD HLLE Riemann Solver
using EnzoRiemannHLLE2 =
  EnzoRiemannImplNew<HLLKernel<EinfeldtWavespeed<HydroLUT>>>;

/// @class    EnzoRiemannHLLMHD2
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLL approximate Riemann Solver using
///           the Davis MHD wavespeed estimator. The same wavespeed estimator
///           is used in Enzo's Runge-Kutta MHD HLL Riemann Solver
using EnzoRiemannHLLMHD2 =EnzoRiemannImplNew<HLLKernel<DavisWavespeed<MHDLUT>>>;

// To define a Riemann Solver for cosmic ray transport, define a separate
// WaveSpeedFunctor struct and define a separate HLL Riemann Solver (maybe
// EnzoRiemannCRHLL)

#endif /* ENZO_ENZO_RIEMANN_HLL2_HPP */
