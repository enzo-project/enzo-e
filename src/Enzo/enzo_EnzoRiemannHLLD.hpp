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
  static int scratch_space_length(const int n_cons_keys) throw(){
    return 2*n_cons_keys;
  }

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
    // This method makes use of the member variables Us and Uss
    // Note that ETA_TOLERANCE is bigger than the tolerance was for the
    // original implementation.

    enzo_float* Us = scratch_space;
    enzo_float* Uss = scratch_space+n_keys;

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
    enzo_float specific_eint_l, specific_eint_r;
    // treating eint as a passively advected scalar for now. In other words,
    //   - specific_eint_l = specific_eint_ls = specific_eint_lss
    //   - specific_eint_r = specific_eint_rs = specific_eint_rss

    enzo_float S_l, S_r, S_ls, S_rs, S_M; // wave speeds
    enzo_float cf_l, cf_r, sam, sap; // fast speeds

    // First, compute Fl and Ul
    rho_l  = prim_l[lut.density];
    p_l    = pressure_l;
    vx_l   = prim_l[lut.velocity_i];
    vy_l   = prim_l[lut.velocity_j];
    vz_l   = prim_l[lut.velocity_k];
    Bx_l   = prim_l[lut.bfield_i];
    By_l   = prim_l[lut.bfield_j];
    Bz_l   = prim_l[lut.bfield_k];
    if (dual_energy){
      specific_eint_l = prim_l[lut.internal_energy];
    }

    Bv_l = Bx_l * vx_l + By_l * vy_l + Bz_l * vz_l;
    etot_l = cons_l[lut.total_energy];
    pt_l = p_l + EnzoRiemann::mag_pressure_(prim_l, lut);
    cf_l = EnzoRiemann::fast_magnetosonic_speed_(prim_l, lut, pressure_l,
						 gamma);

    // load wr and compute the fast magnetosonic speed
    rho_r   = prim_r[lut.density];
    p_r     = pressure_r;
    vx_r    = prim_r[lut.velocity_i];
    vy_r    = prim_r[lut.velocity_j];
    vz_r    = prim_r[lut.velocity_k];
    Bx_r    = prim_r[lut.bfield_i];
    By_r    = prim_r[lut.bfield_j];
    Bz_r    = prim_r[lut.bfield_k];
    if (dual_energy){
      specific_eint_r = prim_r[lut.internal_energy];
    }

    Bv_r = Bx_r * vx_r + By_r * vy_r + Bz_r * vz_r;
    etot_r = cons_r[lut.total_energy];
    pt_r = p_r + EnzoRiemann::mag_pressure_(prim_r, lut);
    cf_r = EnzoRiemann::fast_magnetosonic_speed_(prim_r, lut, pressure_r,
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
      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = flux_l[field];
      }
      vi_bar =  prim_l[lut.velocity_i];
      return;
    } else if (S_r < 0) {
      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = flux_r[field];
      }
      vi_bar =  prim_r[lut.velocity_i];
      return;
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
      setup_cons_ast_(Us, lut, S_M, rho_ls, vy_ls, vz_ls, etot_ls, Bx,
		      By_ls, Bz_ls);

      if (dual_energy){ Us[lut.internal_energy] = rho_ls * specific_eint_l; }

      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = \
	  flux_l[field] + S_l*(Us[field] - cons_l[field]);
      }
      return;
    } else if (S_rs <= 0 && S_r >= 0) {
      // USE F_rs
      setup_cons_ast_(Us, lut, S_M, rho_rs, vy_rs, vz_rs, etot_rs, Bx,
		      By_rs, Bz_rs);

      if (dual_energy){ Us[lut.internal_energy] = rho_rs * specific_eint_r; }

      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = \
	  flux_r[field] + S_r*(Us[field] - cons_r[field]);
      }
      return;
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
      setup_cons_ast_(Us, lut, S_M, rho_ls, vy_ls, vz_ls, etot_ls, Bx,
		      By_ls, Bz_ls);
      setup_cons_ast_(Uss, lut, S_M, rho_ls, vy_ss, vz_ss, etot_lss, Bx,
		      By_ss, Bz_ss);

      if (dual_energy){
	Us[lut.internal_energy]  = rho_ls * specific_eint_l;
	Uss[lut.internal_energy] = rho_ls * specific_eint_l;
      }

      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = \
	  (flux_l[field] + S_ls*Uss[field] - (S_ls - S_l)*Us[field] -
	   S_l*cons_l[field]);
      }
      return;
    } else if (S_M <= 0 && S_rs >= 0) {
      // USE F_rss
      setup_cons_ast_(Us, lut, S_M, rho_rs, vy_rs, vz_rs, etot_rs, Bx,
		      By_rs, Bz_rs);
      setup_cons_ast_(Uss, lut, S_M, rho_rs, vy_ss, vz_ss, etot_rss, Bx,
		      By_ss, Bz_ss);

      if (dual_energy){
	Us[lut.internal_energy]  = rho_rs * specific_eint_r;
	Uss[lut.internal_energy] = rho_rs * specific_eint_r;
      }

      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = \
	  (flux_r[field] + S_rs*Uss[field] - (S_rs - S_r)*Us[field] -
	   S_r*cons_r[field]);
      }
      return;
    }
  }

  static void setup_cons_ast_(enzo_float cons[],
			      const EnzoAdvectionFieldLUT lut,
			      const enzo_float speed, const enzo_float rho,
			      const enzo_float vy, const enzo_float vz,
			      const enzo_float etot, const enzo_float Bx,
			      const enzo_float By, const enzo_float Bz) throw()
  {

    // Helper function that factors out the filling of the of the asterisked
    // and double asterisked conserved quantities
    // The dual energy formalism is explicitly handled outside of this function
    cons[lut.density] = rho;
    cons[lut.velocity_i] = rho * speed;
    cons[lut.velocity_j] = rho * vy;
    cons[lut.velocity_k] = rho * vz;
    cons[lut.total_energy] = etot;
    cons[lut.bfield_i] = Bx;
    cons[lut.bfield_j] = By;
    cons[lut.bfield_k] = Bz;
  }

};


/// @class    EnzoRiemannHLLD
/// @ingroup  Enzo
/// @brief    [\ref Enzo] Encapsulates HLLD approximate Riemann Solver
using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;


#endif /* ENZO_ENZO_RIEMANN_HLLD_HPP */
