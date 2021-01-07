// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannHLLD.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Enzo's HLLD approximate Riemann Solver. Ported from
/// the original Enzo's Riemann_HLLD_MHD.C, written by J. S. Oishi

// this implementation returns identical values to the implementation used in
// Athena for most cases. In the other cases, there are small floating point
// errors.

// Currently, eint fluxes are computed by assuming that specific internal
// energy is a passive scalar. It may be worth considering the calculation of
// fluxes as though it's an actively advected quantity.

#ifndef ENZO_ENZO_RIEMANN_HLLD_HPP
#define ENZO_ENZO_RIEMANN_HLLD_HPP

struct HLLDImpl
{
  /// @class    EnzoRiemannHLLD
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates operations of HLLD approximate Riemann
  /// Solver
public:

  using LUT = EnzoRiemannLUT<MHDLUT>;

  lutarray<LUT> operator()
  (const lutarray<LUT> flux_l, const lutarray<LUT> flux_r,
   const lutarray<LUT> prim_l, const lutarray<LUT> prim_r,
   const lutarray<LUT> cons_l, const lutarray<LUT> cons_r,
   enzo_float pressure_l, enzo_float pressure_r,
   bool barotropic_eos, enzo_float gamma, enzo_float isothermal_cs,
   enzo_float &vi_bar) const noexcept
  {
    // This method makes use of the member variables Us and Uss
    // Note that ETA_TOLERANCE is bigger than the tolerance was for the
    // original implementation.

    lutarray<LUT> Us, Uss, fluxes;

    enzo_float etot_l,etot_r, rho_l, rho_r;
    enzo_float vx_l, vy_l, vz_l, vx_r, vy_r, vz_r;
    enzo_float Bx_l, Bx_r,Bx, By_l, Bz_l, By_r, Bz_r, Bv_l, Bv_r, p_l, p_r;
    enzo_float pt_l, pt_r;
    enzo_float l_coef, r_coef;
    enzo_float rho_ls, rho_rs, vy_ls, vy_rs, vz_ls, vz_rs, vv_ls, vv_rs;
    enzo_float By_ls, By_rs, Bz_ls, Bz_rs, Bv_ls, Bv_rs, bb_ls, bb_rs;
    enzo_float etot_ls, etot_rs, pt_s;
    enzo_float vy_ss, vz_ss, By_ss, Bz_ss, Bv_ss;
    enzo_float etot_lss, etot_rss, rho_savg;

    // when the dual energy formalism is in use, the internal energy is treated
    // as a passively advected scalar and handled separately.
    //
    // In reality, this is not a fully self-consistent treatment of internal
    // energy. Unlike the HLL/HLLC solvers which assumes that the intermediate
    // states have a single constant thermal pressure, the HLLD solver assumes
    // that these states share a single constant total pressure (thermal +
    // magnetic). Thus a more nuanced approach can be taken.
    //
    // Such an approach should be pursued by defining a new LUT class that
    // includes internal energy as an actively advected quantity.

    enzo_float S_l, S_r, S_ls, S_rs, S_M; // wave speeds
    enzo_float cf_l, cf_r, sam, sap; // fast speeds

    // First, compute Fl and Ul
    rho_l  = prim_l[LUT::density];
    p_l    = pressure_l;
    vx_l   = prim_l[LUT::velocity_i];
    vy_l   = prim_l[LUT::velocity_j];
    vz_l   = prim_l[LUT::velocity_k];
    Bx_l   = prim_l[LUT::bfield_i];
    By_l   = prim_l[LUT::bfield_j];
    Bz_l   = prim_l[LUT::bfield_k];

    Bv_l = Bx_l * vx_l + By_l * vy_l + Bz_l * vz_l;
    etot_l = cons_l[LUT::total_energy];
    pt_l = p_l + enzo_riemann_utils::mag_pressure<LUT>(prim_l);
    cf_l = enzo_riemann_utils::fast_magnetosonic_speed<LUT>(prim_l, pressure_l,
                                                            gamma);

    // load wr and compute the fast magnetosonic speed
    rho_r   = prim_r[LUT::density];
    p_r     = pressure_r;
    vx_r    = prim_r[LUT::velocity_i];
    vy_r    = prim_r[LUT::velocity_j];
    vz_r    = prim_r[LUT::velocity_k];
    Bx_r    = prim_r[LUT::bfield_i];
    By_r    = prim_r[LUT::bfield_j];
    Bz_r    = prim_r[LUT::bfield_k];

    Bv_r = Bx_r * vx_r + By_r * vy_r + Bz_r * vz_r;
    etot_r = cons_r[LUT::total_energy];
    pt_r = p_r + enzo_riemann_utils::mag_pressure<LUT>(prim_r);
    cf_r = enzo_riemann_utils::fast_magnetosonic_speed<LUT>(prim_r, pressure_r,
                                                            gamma);

    //
    //wave speeds
    //

    if (Bx_l == Bx_r){
      Bx = Bx_l;
    } else {
      Bx = 0.5*(Bx_l + Bx_r);
    }
    // first, outermost wave speeds
    // simplest choice from Miyoshi & Kusano (2005)
    S_l = std::min(vx_l, vx_r) - std::max(cf_l, cf_r);
    S_r = std::max(vx_l, vx_r) + std::max(cf_l, cf_r);

    // if S_l>0 or S_r<0, no need to go further (the internal energy is also
    // automatically handled if the dual energy formalism is in use)
    if (S_l > 0) {
      for (std::size_t field = 0; field<LUT::NEQ; field++){
	fluxes[field] = flux_l[field];
      }
      vi_bar =  prim_l[LUT::velocity_i];
      return fluxes;
    } else if (S_r < 0) {
      for (std::size_t field = 0; field<LUT::NEQ; field++){
	fluxes[field] = flux_r[field];
      }
      vi_bar =  prim_r[LUT::velocity_i];
      return fluxes;
    } 

    // next, the middle (contact) wave
    S_M = ((S_r - vx_r)*rho_r*vx_r - (S_l - vx_l)*
	   rho_l*vx_l - pt_r + pt_l)/((S_r - vx_r)*
				      rho_r - (S_l - vx_l)*rho_l);

    // Brief segue to compute vx (velocity component normal to the interface).
    // Following the convention from Enzo's flux_hllc.F:
    //     - when S_l > 0 use vi_bar = vx_L [see above]
    //     - when S_r < 0 use vi_bar = vx_R [see above]
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
    l_coef = (S_l - vx_l)/(S_l - S_M);
    r_coef = (S_r - vx_r)/(S_r - S_M);
    if (S_M >=0){
      vi_bar = S_M * l_coef;
    } else {
      vi_bar = S_M * r_coef;
    }

    // finally, compute the sppeds of the intermediate (Alfven) waves
    rho_ls = rho_l * l_coef;
    rho_rs = rho_r * r_coef;

    S_ls = S_M - std::fabs(Bx)/std::sqrt(rho_ls);
    S_rs = S_M + std::fabs(Bx)/std::sqrt(rho_rs);

    pt_s = ((S_r -  vx_r) * rho_r*pt_l - (S_l - vx_l) * rho_l * pt_r +
	    rho_l*rho_r*(S_r - vx_r)*(S_l - vx_l)*
	    (vx_r - vx_l))/((S_r - vx_r)*rho_r - (S_l - vx_l)*rho_l);

    sam = vx_l - cf_l;
    sap = vx_l + cf_l;
      
    if ((std::fabs(S_M - vx_l) <= ETA_TOLERANCE) and 
        (std::fabs(By_l) <= ETA_TOLERANCE) and 
        (std::fabs(Bz_l) <= ETA_TOLERANCE) and 
        (Bx*Bx >= gamma * p_l) and
        ((std::fabs(S_l - sam) <= ETA_TOLERANCE) or
	 (std::fabs(S_l - sap) <= ETA_TOLERANCE)) ) {
      vy_ls = vy_l;
      vz_ls = vz_l;
      By_ls = By_l;
      Bz_ls = Bz_l;
    } else {
      vv_ls = (S_M - vx_l)/(rho_l*(S_l - vx_l)*(S_l - S_M) - Bx*Bx);
      bb_ls = ((rho_l*(S_l - vx_l)*(S_l - vx_l) - Bx*Bx) /
	       (rho_l*(S_l - vx_l)*(S_l - S_M) - Bx*Bx));
      vy_ls = vy_l - Bx * By_l * vv_ls;
      By_ls = By_l * bb_ls;
      vz_ls = vz_l - Bx * Bz_l * vv_ls;
      Bz_ls = Bz_l * bb_ls;
    }

    sam = vx_r - cf_r;
    sap = vx_r + cf_r;

    if ((std::fabs(S_M - vx_r) <= ETA_TOLERANCE) and 
        (std::fabs(By_r) <= ETA_TOLERANCE) and 
        (std::fabs(Bz_r) <= ETA_TOLERANCE) and 
        (Bx*Bx >= gamma * p_r) and
        ((std::fabs(S_r - sam) <= ETA_TOLERANCE) or
	 (std::fabs(S_r - sap) <= ETA_TOLERANCE)) ) {
      vy_rs = vy_r;
      vz_rs = vz_r;
      By_rs = By_r;
      Bz_rs = Bz_r;
    } else {
      vv_rs = (S_M - vx_r)/(rho_r*(S_r - vx_r)*(S_r - S_M) - Bx*Bx);
      bb_rs = ((rho_r*(S_r - vx_r)*(S_r - vx_r) - Bx*Bx)/
	       (rho_r*(S_r - vx_r)*(S_r - S_M) - Bx*Bx));
      vy_rs = vy_r - Bx * By_r * vv_rs;
      vz_rs = vz_r - Bx * Bz_r * vv_rs;
      By_rs = By_r * bb_rs;
      Bz_rs = Bz_r * bb_rs;
    }
    Bv_ls = S_M * Bx + vy_ls * By_ls + vz_ls * Bz_ls;
    Bv_rs = S_M * Bx + vy_rs * By_rs + vz_rs * Bz_rs;

    etot_ls = ((S_l - vx_l)*etot_l - pt_l*vx_l + pt_s * S_M +
	       Bx*(Bv_l - Bv_ls))/(S_l - S_M);
    etot_rs = ((S_r - vx_r)*etot_r - pt_r*vx_r + pt_s * S_M +
	       Bx*(Bv_r - Bv_rs))/(S_r - S_M);

    // compute the fluxes based on the wave speeds
    if (S_l <= 0 && S_ls >= 0) {
      // USE F_ls
      Us = setup_cons_ast_(S_M, rho_ls, vy_ls, vz_ls, etot_ls, Bx,
			   By_ls, Bz_ls);

      for (std::size_t field = 0; field<LUT::NEQ; field++){
	fluxes[field] = flux_l[field] + S_l*(Us[field] - cons_l[field]);
      }
      return fluxes;
    } else if (S_rs <= 0 && S_r >= 0) {
      // USE F_rs
      Us = setup_cons_ast_(S_M, rho_rs, vy_rs, vz_rs, etot_rs, Bx,
			   By_rs, Bz_rs);

      for (std::size_t field = 0; field<LUT::NEQ; field++){
	fluxes[field] = flux_r[field] + S_r*(Us[field] - cons_r[field]);
      }
      return fluxes;
    }

    enzo_float sign_Bx = (Bx>0) ? 1. : -1.;

    //do U** stuff
    rho_savg = std::sqrt(rho_ls) + std::sqrt(rho_rs);
    vy_ss = (std::sqrt(rho_ls) * vy_ls + std::sqrt(rho_rs) * vy_rs +
	     (By_rs - By_ls) * sign_Bx)/rho_savg;
    vz_ss = (std::sqrt(rho_ls) * vz_ls + std::sqrt(rho_rs) * vz_rs +
	     (Bz_rs - Bz_ls) * sign_Bx)/rho_savg;
    By_ss = (std::sqrt(rho_ls) * By_rs + std::sqrt(rho_rs) * By_ls +
	     std::sqrt(rho_ls * rho_rs) * (vy_rs - vy_ls) * sign_Bx)/rho_savg;
    Bz_ss = (std::sqrt(rho_ls) * Bz_rs + std::sqrt(rho_rs) * Bz_ls +
	     std::sqrt(rho_ls * rho_rs) * (vz_rs - vz_ls) * sign_Bx)/rho_savg;
    Bv_ss = S_M * Bx + vy_ss * By_ss + vz_ss * Bz_ss;
    etot_lss = etot_ls - std::sqrt(rho_ls) * (Bv_ls - Bv_ss) * sign_Bx;
    etot_rss = etot_rs + std::sqrt(rho_rs) * (Bv_rs - Bv_ss) * sign_Bx;

    if (S_ls <= 0 && S_M >= 0) {
      // USE F_lss
      Us = setup_cons_ast_(S_M, rho_ls, vy_ls, vz_ls, etot_ls,
			   Bx, By_ls, Bz_ls);
      Uss = setup_cons_ast_(S_M, rho_ls, vy_ss, vz_ss, etot_lss,
			    Bx, By_ss, Bz_ss);

      for (std::size_t field = 0; field<LUT::NEQ; field++){
	fluxes[field] = (flux_l[field] + S_ls*Uss[field] -
			 (S_ls - S_l)*Us[field] - S_l*cons_l[field]);
      }
      return fluxes;
    } else if (S_M <= 0 && S_rs >= 0) {
      // USE F_rss
      Us = setup_cons_ast_(S_M, rho_rs, vy_rs, vz_rs, etot_rs,
			   Bx, By_rs, Bz_rs);
      Uss = setup_cons_ast_(S_M, rho_rs, vy_ss, vz_ss, etot_rss,
			    Bx, By_ss, Bz_ss);

      for (std::size_t field = 0; field<LUT::NEQ; field++){
	fluxes[field] = (flux_r[field] + S_rs*Uss[field] -
			 (S_rs - S_r)*Us[field] - S_r*cons_r[field]);
      }
      return fluxes;
    }
    // include this to avoid getting yelled at
    return fluxes;
  }

  lutarray<LUT> setup_cons_ast_( const enzo_float speed,
			       const enzo_float rho, const enzo_float vy,
			       const enzo_float vz, const enzo_float etot,
			       const enzo_float Bx, const enzo_float By,
			       const enzo_float Bz) const noexcept
  {

    // Helper function that factors out the filling of the of the asterisked
    // and double asterisked conserved quantities
    // The dual energy formalism is explicitly handled outside of this function
    lutarray<LUT> out;
    out[LUT::density] = rho;
    out[LUT::velocity_i] = rho * speed;
    out[LUT::velocity_j] = rho * vy;
    out[LUT::velocity_k] = rho * vz;
    out[LUT::total_energy] = etot;
    out[LUT::bfield_i] = Bx;
    out[LUT::bfield_j] = By;
    out[LUT::bfield_k] = Bz;
    return out;
  }

};


/// @class    EnzoRiemannHLLD
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLD approximate Riemann Solver
using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;


#endif /* ENZO_ENZO_RIEMANN_HLLD_HPP */
