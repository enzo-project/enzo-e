// See LICENSE_ATHENAPP file for license and copyright information

/// @file     enzo_EnzoRiemannHLLD.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Enzo's HLLD approximate Riemann Solver. Ported from
/// Athena++

// Currently, eint fluxes are computed by assuming that specific internal
// energy is a passive scalar. It may be worth considering the calculation of
// fluxes as though it's an actively advected quantity.

#ifndef ENZO_ENZO_RIEMANN_HLLD_HPP
#define ENZO_ENZO_RIEMANN_HLLD_HPP

#ifndef SQR
#define SQR(x) ( (x)*(x) )
#endif

#define SMALL_NUMBER 1.0e-8

struct HLLDKernel
{
  /// @class    HLLDKernel
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of HLLD approximate Riemann
  /// Solver

public: // typedefs
  using LUT = EnzoRiemannLUT<MHDLUT>;
  using EOSStructT = EOSStructIdeal;

  struct Cons1D { enzo_float d, mx, my, mz, e, by, bz; };

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
    const int external_bfield_i = config.dim + LUT::bfield_i;
    const int external_bfield_j = ((config.dim+1)%3) + LUT::bfield_i;
    const int external_bfield_k = ((config.dim+2)%3) + LUT::bfield_i;

    const enzo_float gamma = config.eos.get_gamma();
    const enzo_float igm1 = 1.0 / (gamma - 1.0);

    lutarray<LUT> flxi;      // temporary variable to store flux
    lutarray<LUT> wli, wri;  // Left, Right reconstructed primitive
    enzo_float spd[5];       // signal speeds, left to right

    Cons1D ul,ur;                 // L/R states, conserved variables (computed)
    Cons1D ulst,uldst,urdst,urst; // Conserved variable for all states
    Cons1D fl,fr;                 // Fluxes for left & right states

    //--- Step 1.  Load L/R states into local variables

    wli[LUT::density]      = config.prim_arr_l(LUT::density,iz,iy,ix);
    wli[LUT::velocity_i]   = config.prim_arr_l(external_velocity_i,iz,iy,ix);
    wli[LUT::velocity_j]   = config.prim_arr_l(external_velocity_j,iz,iy,ix);
    wli[LUT::velocity_k]   = config.prim_arr_l(external_velocity_k,iz,iy,ix);
    wli[LUT::bfield_i]     = config.prim_arr_l(external_bfield_i,iz,iy,ix);
    wli[LUT::bfield_j]     = config.prim_arr_l(external_bfield_j,iz,iy,ix);
    wli[LUT::bfield_k]     = config.prim_arr_l(external_bfield_k,iz,iy,ix);
    // the following actually stores pressure, not total_energy
    wli[LUT::total_energy] = config.prim_arr_l(LUT::total_energy,iz,iy,ix);

    wri[LUT::density]      = config.prim_arr_r(LUT::density,iz,iy,ix);
    wri[LUT::velocity_i]   = config.prim_arr_r(external_velocity_i,iz,iy,ix);
    wri[LUT::velocity_j]   = config.prim_arr_r(external_velocity_j,iz,iy,ix);
    wri[LUT::velocity_k]   = config.prim_arr_r(external_velocity_k,iz,iy,ix);
    wri[LUT::bfield_i]     = config.prim_arr_r(external_bfield_i,iz,iy,ix);
    wri[LUT::bfield_j]     = config.prim_arr_r(external_bfield_j,iz,iy,ix);
    wri[LUT::bfield_k]     = config.prim_arr_r(external_bfield_k,iz,iy,ix);
    // the following actually stores pressure, not total_energy
    wri[LUT::total_energy] = config.prim_arr_r(LUT::total_energy,iz,iy,ix);


    // load left and right pressure values
    const enzo_float pressure_l = wli[LUT::total_energy];
    const enzo_float pressure_r = wri[LUT::total_energy];

    const enzo_float bxi = wli[LUT::bfield_i]; // == wri[LUT::bfield_i]

    // Compute L/R states for selected conserved variables
    enzo_float bxsq = bxi*bxi;
    // group transverse vector components for floating-point associativity
    // symmetry
    // magnetic pressure (l/r):
    enzo_float pbl = 0.5*(bxsq + (SQR(wli[LUT::bfield_j]) +
                                  SQR(wli[LUT::bfield_k])));
    enzo_float pbr = 0.5*(bxsq + (SQR(wri[LUT::bfield_j]) +
                                  SQR(wri[LUT::bfield_k])));
    // kinetic energy:
    enzo_float kel = 0.5*wli[LUT::density]*(SQR(wli[LUT::velocity_i]) +
                                            (SQR(wli[LUT::velocity_j]) +
                                             SQR(wli[LUT::velocity_k])));
    enzo_float ker = 0.5*wri[LUT::density]*(SQR(wri[LUT::velocity_i]) +
                                            (SQR(wri[LUT::velocity_j]) +
                                             SQR(wri[LUT::velocity_k])));

    ul.d  = wli[LUT::density];
    ul.mx = wli[LUT::velocity_i]*ul.d;
    ul.my = wli[LUT::velocity_j]*ul.d;
    ul.mz = wli[LUT::velocity_k]*ul.d;
    ul.e  = pressure_l*igm1 + kel + pbl;
    ul.by = wli[LUT::bfield_j];
    ul.bz = wli[LUT::bfield_k];

    ur.d  = wri[LUT::density];
    ur.mx = wri[LUT::velocity_i]*ur.d;
    ur.my = wri[LUT::velocity_j]*ur.d;
    ur.mz = wri[LUT::velocity_k]*ur.d;
    ur.e  = pressure_r*igm1 + ker + pbr;
    ur.by = wri[LUT::bfield_j];
    ur.bz = wri[LUT::bfield_k];

    //--- Step 2.  Compute L & R wave speeds according to Miyoshi & Kusano, eqn. (67)

    using enzo_riemann_utils::fast_magnetosonic_speed;
    enzo_float cfl = fast_magnetosonic_speed<LUT>(wli, pressure_l, config.eos);
    enzo_float cfr = fast_magnetosonic_speed<LUT>(wri, pressure_r, config.eos);

    spd[0] = std::min( wli[LUT::velocity_i]-cfl, wri[LUT::velocity_i]-cfr );
    spd[4] = std::max( wli[LUT::velocity_i]+cfl, wri[LUT::velocity_i]+cfr );

    // Code analogous to the following block was commented out:
    // enzo_float cfmax = std::max(cfl,cfr);
    // if (wli[LUT::velocity_i] <= wri[LUT::velocity_i]) {
    //   spd[0] = wli[LUT::velocity_i] - cfmax;
    //   spd[4] = wri[LUT::velocity_i] + cfmax;
    // } else {
    //   spd[0] = wri[LUT::velocity_i] - cfmax;
    //   spd[4] = wli[LUT::velocity_i] + cfmax;
    // }

    //--- Step 3.  Compute L/R fluxes

    enzo_float ptl = pressure_l + pbl; // total pressures L,R
    enzo_float ptr = pressure_r + pbr;

    fl.d  = ul.mx;
    fl.mx = ul.mx*wli[LUT::velocity_i] + ptl - bxsq;
    fl.my = ul.my*wli[LUT::velocity_i] - bxi*ul.by;
    fl.mz = ul.mz*wli[LUT::velocity_i] - bxi*ul.bz;
    fl.e  = wli[LUT::velocity_i]*(ul.e + ptl - bxsq) - bxi*(wli[LUT::velocity_j]*ul.by +
                                                            wli[LUT::velocity_k]*ul.bz);
    fl.by = ul.by*wli[LUT::velocity_i] - bxi*wli[LUT::velocity_j];
    fl.bz = ul.bz*wli[LUT::velocity_i] - bxi*wli[LUT::velocity_k];

    fr.d  = ur.mx;
    fr.mx = ur.mx*wri[LUT::velocity_i] + ptr - bxsq;
    fr.my = ur.my*wri[LUT::velocity_i] - bxi*ur.by;
    fr.mz = ur.mz*wri[LUT::velocity_i] - bxi*ur.bz;
    fr.e  = wri[LUT::velocity_i]*(ur.e + ptr - bxsq) - bxi*(wri[LUT::velocity_j]*ur.by +
                                                            wri[LUT::velocity_k]*ur.bz);
    fr.by = ur.by*wri[LUT::velocity_i] - bxi*wri[LUT::velocity_j];
    fr.bz = ur.bz*wri[LUT::velocity_i] - bxi*wri[LUT::velocity_k];

    //--- Step 4.  Compute middle and Alfven wave speeds

    enzo_float sdl = spd[0] - wli[LUT::velocity_i];  // S_i-u_i (i=L or R)
    enzo_float sdr = spd[4] - wri[LUT::velocity_i];

    // S_M: eqn (38) of Miyoshi & Kusano
    // (KGF): group ptl, ptr terms for floating-point associativity symmetry
    spd[2] = (sdr*ur.mx - sdl*ul.mx + (ptl - ptr))/(sdr*ur.d - sdl*ul.d);

    enzo_float sdml   = spd[0] - spd[2];  // S_i-S_M (i=L or R)
    enzo_float sdmr   = spd[4] - spd[2];
    enzo_float sdml_inv = 1.0/sdml;
    enzo_float sdmr_inv = 1.0/sdmr;
    // eqn (43) of Miyoshi & Kusano
    ulst.d = ul.d * sdl * sdml_inv;
    urst.d = ur.d * sdr * sdmr_inv;
    enzo_float ulst_d_inv = 1.0/ulst.d;
    enzo_float urst_d_inv = 1.0/urst.d;
    enzo_float sqrtdl = std::sqrt(ulst.d);
    enzo_float sqrtdr = std::sqrt(urst.d);

    // eqn (51) of Miyoshi & Kusano
    spd[1] = spd[2] - std::abs(bxi)/sqrtdl;
    spd[3] = spd[2] + std::abs(bxi)/sqrtdr;

    //--- Step 5.  Compute intermediate states
    // eqn (23) explicitly becomes eq (41) of Miyoshi & Kusano
    // place an assertion that ptstl==ptstr
    enzo_float ptstl = ptl + ul.d*sdl*(spd[2]-wli[LUT::velocity_i]);
    enzo_float ptstr = ptr + ur.d*sdr*(spd[2]-wri[LUT::velocity_i]);
    // the following 2 equations had issues when averaged
    // enzo_float ptstl = ptl + ul.d*sdl*(sdl-sdml); 
    // enzo_float ptstr = ptr + ur.d*sdr*(sdr-sdmr);
    enzo_float ptst = 0.5*(ptstr + ptstl);  // total pressure (star state)

    // ul* - eqn (39) of M&K
    ulst.mx = ulst.d * spd[2];
    if (std::abs(ul.d*sdl*sdml-bxsq) < (SMALL_NUMBER)*ptst) {
      // Degenerate case
      ulst.my = ulst.d * wli[LUT::velocity_j];
      ulst.mz = ulst.d * wli[LUT::velocity_k];

      ulst.by = ul.by;
      ulst.bz = ul.bz;
    } else {
      // eqns (44) and (46) of M&K
      enzo_float tmp = bxi*(sdl - sdml)/(ul.d*sdl*sdml - bxsq);
      ulst.my = ulst.d * (wli[LUT::velocity_j] - ul.by*tmp);
      ulst.mz = ulst.d * (wli[LUT::velocity_k] - ul.bz*tmp);

      // eqns (45) and (47) of M&K
      tmp = (ul.d*SQR(sdl) - bxsq)/(ul.d*sdl*sdml - bxsq);
      ulst.by = ul.by * tmp;
      ulst.bz = ul.bz * tmp;
    }
    // v_i* dot B_i*
    // (KGF): group transverse momenta terms for floating-point associativity symmetry
    enzo_float vbstl = (ulst.mx*bxi + (ulst.my*ulst.by +
                                       ulst.mz*ulst.bz)) * ulst_d_inv;
    // eqn (48) of M&K
    // (KGF): group transverse by, bz terms for floating-point associativity symmetry
    ulst.e = (sdl*ul.e - ptl*wli[LUT::velocity_i] + ptst*spd[2] +
              bxi*(wli[LUT::velocity_i]*bxi +
                   (wli[LUT::velocity_j]*ul.by + wli[LUT::velocity_k]*ul.bz)
                   - vbstl))*sdml_inv;

    // ur* - eqn (39) of M&K
    urst.mx = urst.d * spd[2];
    if (std::abs(ur.d*sdr*sdmr - bxsq) < (SMALL_NUMBER)*ptst) {
      // Degenerate case
      urst.my = urst.d * wri[LUT::velocity_j];
      urst.mz = urst.d * wri[LUT::velocity_k];

      urst.by = ur.by;
      urst.bz = ur.bz;
    } else {
      // eqns (44) and (46) of M&K
      enzo_float tmp = bxi*(sdr - sdmr)/(ur.d*sdr*sdmr - bxsq);
      urst.my = urst.d * (wri[LUT::velocity_j] - ur.by*tmp);
      urst.mz = urst.d * (wri[LUT::velocity_k] - ur.bz*tmp);

      // eqns (45) and (47) of M&K
      tmp = (ur.d*SQR(sdr) - bxsq)/(ur.d*sdr*sdmr - bxsq);
      urst.by = ur.by * tmp;
      urst.bz = ur.bz * tmp;
    }
    // v_i* dot B_i*
    // (KGF): group transverse momenta terms for floating-point associativity symmetry
    enzo_float vbstr = (urst.mx*bxi+(urst.my*urst.by+
                                     urst.mz*urst.bz))*urst_d_inv;
    // eqn (48) of M&K
    // (KGF): group transverse by, bz terms for floating-point associativity symmetry
    urst.e = (sdr*ur.e - ptr*wri[LUT::velocity_i] + ptst*spd[2] +
              bxi*(wri[LUT::velocity_i]*bxi +
                   (wri[LUT::velocity_j]*ur.by +
                    wri[LUT::velocity_k]*ur.bz) - vbstr))*sdmr_inv;
    // ul** and ur** - if Bx is near zero, same as *-states
    if (0.5*bxsq < (SMALL_NUMBER)*ptst) {
      uldst = ulst;
      urdst = urst;
    } else {
      enzo_float invsumd = 1.0/(sqrtdl + sqrtdr);
      enzo_float bxsig = (bxi > 0.0 ? 1.0 : -1.0);

      uldst.d = ulst.d;
      urdst.d = urst.d;

      uldst.mx = ulst.mx;
      urdst.mx = urst.mx;

      // eqn (59) of M&K
      enzo_float tmp = invsumd*(sqrtdl*(ulst.my*ulst_d_inv) +
                                sqrtdr*(urst.my*urst_d_inv) +
                                bxsig*(urst.by - ulst.by));
      uldst.my = uldst.d * tmp;
      urdst.my = urdst.d * tmp;

      // eqn (60) of M&K
      tmp = invsumd*(sqrtdl*(ulst.mz*ulst_d_inv) +
                     sqrtdr*(urst.mz*urst_d_inv) +
                     bxsig*(urst.bz - ulst.bz));
      uldst.mz = uldst.d * tmp;
      urdst.mz = urdst.d * tmp;

      // eqn (61) of M&K
      tmp = invsumd*(sqrtdl*urst.by + sqrtdr*ulst.by +
                     bxsig*sqrtdl*sqrtdr*((urst.my*urst_d_inv) -
                                          (ulst.my*ulst_d_inv)));
      uldst.by = urdst.by = tmp;

      // eqn (62) of M&K
      tmp = invsumd*(sqrtdl*urst.bz + sqrtdr*ulst.bz +
                     bxsig*sqrtdl*sqrtdr*((urst.mz*urst_d_inv) -
                                          (ulst.mz*ulst_d_inv)));
      uldst.bz = urdst.bz = tmp;

      // eqn (63) of M&K
      tmp = spd[2]*bxi + (uldst.my*uldst.by + uldst.mz*uldst.bz)/uldst.d;
      uldst.e = ulst.e - sqrtdl*bxsig*(vbstl - tmp);
      urdst.e = urst.e + sqrtdr*bxsig*(vbstr - tmp);
    }

    //--- Step 6.  Compute flux
    uldst.d = spd[1] * (uldst.d - ulst.d);
    uldst.mx = spd[1] * (uldst.mx - ulst.mx);
    uldst.my = spd[1] * (uldst.my - ulst.my);
    uldst.mz = spd[1] * (uldst.mz - ulst.mz);
    uldst.e = spd[1] * (uldst.e - ulst.e);
    uldst.by = spd[1] * (uldst.by - ulst.by);
    uldst.bz = spd[1] * (uldst.bz - ulst.bz);

    ulst.d = spd[0] * (ulst.d - ul.d);
    ulst.mx = spd[0] * (ulst.mx - ul.mx);
    ulst.my = spd[0] * (ulst.my - ul.my);
    ulst.mz = spd[0] * (ulst.mz - ul.mz);
    ulst.e = spd[0] * (ulst.e - ul.e);
    ulst.by = spd[0] * (ulst.by - ul.by);
    ulst.bz = spd[0] * (ulst.bz - ul.bz);

    urdst.d = spd[3] * (urdst.d - urst.d);
    urdst.mx = spd[3] * (urdst.mx - urst.mx);
    urdst.my = spd[3] * (urdst.my - urst.my);
    urdst.mz = spd[3] * (urdst.mz - urst.mz);
    urdst.e = spd[3] * (urdst.e - urst.e);
    urdst.by = spd[3] * (urdst.by - urst.by);
    urdst.bz = spd[3] * (urdst.bz - urst.bz);

    urst.d = spd[4] * (urst.d  - ur.d);
    urst.mx = spd[4] * (urst.mx - ur.mx);
    urst.my = spd[4] * (urst.my - ur.my);
    urst.mz = spd[4] * (urst.mz - ur.mz);
    urst.e = spd[4] * (urst.e - ur.e);
    urst.by = spd[4] * (urst.by - ur.by);
    urst.bz = spd[4] * (urst.bz - ur.bz);

    if (spd[0] >= 0.0) {
      // return Fl if flow is supersonic
      flxi[LUT::density] = fl.d;
      flxi[LUT::velocity_i] = fl.mx;
      flxi[LUT::velocity_j] = fl.my;
      flxi[LUT::velocity_k] = fl.mz;
      flxi[LUT::total_energy] = fl.e;
      flxi[LUT::bfield_j] = fl.by;
      flxi[LUT::bfield_k] = fl.bz;
    } else if (spd[4] <= 0.0) {
      // return Fr if flow is supersonic
      flxi[LUT::density] = fr.d;
      flxi[LUT::velocity_i] = fr.mx;
      flxi[LUT::velocity_j] = fr.my;
      flxi[LUT::velocity_k] = fr.mz;
      flxi[LUT::total_energy] = fr.e;
      flxi[LUT::bfield_j] = fr.by;
      flxi[LUT::bfield_k] = fr.bz;
    } else if (spd[1] >= 0.0) {
      // return Fl*
      flxi[LUT::density] = fl.d  + ulst.d;
      flxi[LUT::velocity_i] = fl.mx + ulst.mx;
      flxi[LUT::velocity_j] = fl.my + ulst.my;
      flxi[LUT::velocity_k] = fl.mz + ulst.mz;
      flxi[LUT::total_energy] = fl.e  + ulst.e;
      flxi[LUT::bfield_j] = fl.by + ulst.by;
      flxi[LUT::bfield_k] = fl.bz + ulst.bz;
    } else if (spd[2] >= 0.0) {
      // return Fl**
      flxi[LUT::density] = fl.d  + ulst.d + uldst.d;
      flxi[LUT::velocity_i] = fl.mx + ulst.mx + uldst.mx;
      flxi[LUT::velocity_j] = fl.my + ulst.my + uldst.my;
      flxi[LUT::velocity_k] = fl.mz + ulst.mz + uldst.mz;
      flxi[LUT::total_energy] = fl.e  + ulst.e + uldst.e;
      flxi[LUT::bfield_j] = fl.by + ulst.by + uldst.by;
      flxi[LUT::bfield_k] = fl.bz + ulst.bz + uldst.bz;
    } else if (spd[3] > 0.0) {
      // return Fr**
      flxi[LUT::density] = fr.d + urst.d + urdst.d;
      flxi[LUT::velocity_i] = fr.mx + urst.mx + urdst.mx;
      flxi[LUT::velocity_j] = fr.my + urst.my + urdst.my;
      flxi[LUT::velocity_k] = fr.mz + urst.mz + urdst.mz;
      flxi[LUT::total_energy] = fr.e + urst.e + urdst.e;
      flxi[LUT::bfield_j] = fr.by + urst.by + urdst.by;
      flxi[LUT::bfield_k] = fr.bz + urst.bz + urdst.bz;
    } else {
      // return Fr*
      flxi[LUT::density] = fr.d  + urst.d;
      flxi[LUT::velocity_i] = fr.mx + urst.mx;
      flxi[LUT::velocity_j] = fr.my + urst.my;
      flxi[LUT::velocity_k] = fr.mz + urst.mz;
      flxi[LUT::total_energy] = fr.e  + urst.e;
      flxi[LUT::bfield_j] = fr.by + urst.by;
      flxi[LUT::bfield_k] = fr.bz + urst.bz;
    }

    config.flux_arr(LUT::density,iz,iy,ix) = flxi[LUT::density];
    config.flux_arr(external_velocity_i,iz,iy,ix) = flxi[LUT::velocity_i];
    config.flux_arr(external_velocity_j,iz,iy,ix) = flxi[LUT::velocity_j];
    config.flux_arr(external_velocity_k,iz,iy,ix) = flxi[LUT::velocity_k];
    config.flux_arr(external_bfield_i,iz,iy,ix) = 0.0;
    config.flux_arr(external_bfield_j,iz,iy,ix) = flxi[LUT::bfield_j];
    config.flux_arr(external_bfield_k,iz,iy,ix) = flxi[LUT::bfield_k];
    config.flux_arr(LUT::total_energy,iz,iy,ix) = flxi[LUT::total_energy];

    // finally, deal with dual energy stuff.
    // compute internal energy flux, assuming passive advection
    // (this was not handled with the rest of the fluxes)
    config.internal_energy_flux_arr(iz,iy,ix) =
      enzo_riemann_utils::passive_eint_flux
      (wli[LUT::density], pressure_l, wri[LUT::density], pressure_r,
       config.eos, config.flux_arr(LUT::density,iz,iy,ix));

    // compute vi_bar, velocity component normal to the interface
    // for simplicity, we adopt the shorthand:
    //     - wli[LUT::velocity_i] <--> vx_L
    //     - wri[LUT::velocity_i] <--> vx_R
    // 
    // Following the convention from Enzo's flux_hllc.F:
    //     - when S_l > 0 use vi_bar = vx_L
    //     - when S_r < 0 use vi_bar = vx_R
    //     - otherwise, linearly interpolate the velocity at x=0 (the cell
    //       interface) at time t (which cancels out) between the (x,v) points:
    //         (S_l*t, vx_L) and (S_M*t, S_M)  if S_l <= 0 & S_M >= 0
    //         (S_M*t, S_M)  and (S_r*t, vx_R) if S_M <= 0 & S_r >= 0
    //       recall that S_l <= S_M <= S_r.
    // Note that the last case isn't completely consistent with the underlying
    // assumption of the HLLD (and HLLC) solver that vx is constant throughout
    // the full intermediate region between S_l and S_r and equal to S_M. To be
    // completely consistent, we should just use S_M. This lack of consistency
    // also explains why we don't worry about the intermediate alfven waves for
    // computing vi_bar (they have definition of vx other than vx = S_M).
    const enzo_float S_M = spd[2];
    const enzo_float S_l = spd[0];
    const enzo_float S_r = spd[4];

    const enzo_float l_coef = (S_l - wli[LUT::velocity_i])/(S_l - S_M);
    const enzo_float r_coef = (S_r - wri[LUT::velocity_i])/(S_r - S_M);
    if (S_l > 0) {
      config.velocity_i_bar_arr(iz,iy,ix) = wli[LUT::velocity_i];
    } else if (S_r < 0) {
      config.velocity_i_bar_arr(iz,iy,ix) = wri[LUT::velocity_i];
    } else if (S_M >=0){
      config.velocity_i_bar_arr(iz,iy,ix) = S_M * l_coef;
    } else {
      config.velocity_i_bar_arr(iz,iy,ix) = S_M * r_coef;
    }
  }
};

/// @class    EnzoRiemannHLLD
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLD approximate Riemann Solver
using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDKernel>;

#endif /* ENZO_ENZO_RIEMANN_HLLD_HPP */
