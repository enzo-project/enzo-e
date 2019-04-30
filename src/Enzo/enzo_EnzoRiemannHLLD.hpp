#ifndef ENZO_ENZO_RIEMANN_HLLD_HPP
#define ENZO_ENZO_RIEMANN_HLLD_HPP

// This is a quick port from Enzo's Riemann_HLLD_MHD.C, written by J. S. Oishi
// The current implementation assumes that each process will always have an
// independent instance of this object. If that changes, then handling of the
// private members Us and Uss will also need to change.

class EnzoRiemannHLLD : public EnzoRiemann
{

  /// @class    EnzoRiemannHLLD
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates HLLE approximate Riemann Solver

public: // interface
  /// Create a new EnzoRiemannHLLD object
  EnzoRiemannHLLD(EnzoFieldConditions cond,
		  std::vector<std::string> &passive_groups,
		  FluxFunctor** flux_funcs, int n_funcs)
    : EnzoRiemann(cond, passive_groups, flux_funcs, n_funcs)
  {
    // Reserve space in the following scratch space arrays
    Us = new enzo_float[this->n_cons_keys_];
    Uss = new enzo_float[this->n_cons_keys_];
  }

  /// Virtual destructor
  virtual ~EnzoRiemannHLLD()
  {
    delete[] Us;
    delete[] Uss;
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoRiemannHLLD);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemannHLLD (CkMigrateMessage *m)
    : EnzoRiemann(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    EnzoRiemann::pup(p);
    // Reserve space in the following scratch space arrays
    Us = new enzo_float[this->n_cons_keys_];
    Uss = new enzo_float[this->n_cons_keys_];
  };

  void calc_riemann_fluxes_(const enzo_float flux_l[],
			    const enzo_float flux_r[],
			    const enzo_float prim_l[],
			    const enzo_float prim_r[],
			    const enzo_float cons_l[],
			    const enzo_float cons_r[],
			    const field_lut prim_lut,
			    const field_lut cons_lut,
			    const int n_keys,
			    EnzoEquationOfState *eos,
			    const int iz, const int iy, const int ix,
			    EFlt3DArray flux_arrays[])
  {
    // This method makes use of the member variables Us and Uss
    // Note that ETA_TOLERANCE is bigger than the tolerance was for the
    // original implementation.

    enzo_float etot_l,etot_r, rho_l, rho_r;
    enzo_float vx_l, vy_l, vz_l, vx_r, vy_r, vz_r;
    enzo_float Bx_l, Bx_r,Bx, By_l, Bz_l, By_r, Bz_r, Bv_l, Bv_r, p_l, p_r;
    enzo_float pt_l, pt_r;
    enzo_float rho_ls, rho_rs, vy_ls, vy_rs, vz_ls, vz_rs, vv_ls, vv_rs;
    enzo_float By_ls, By_rs, Bz_ls, Bz_rs, Bv_ls, Bv_rs, bb_ls, bb_rs;
    enzo_float etot_ls, etot_rs, pt_s;
    enzo_float vy_ss, vz_ss, By_ss, Bz_ss, Bv_ss;
    enzo_float etot_lss, etot_rss, rho_savg;

    enzo_float S_l, S_r, S_ls, S_rs, S_M; // wave speeds
    enzo_float cf_l, cf_r, sam, sap; // fast speeds

    enzo_float gamma = eos->get_gamma();

    // First, compute Fl and Ul
    rho_l  = prim_l[prim_lut.density];
    p_l    = prim_l[prim_lut.pressure];
    vx_l   = prim_l[prim_lut.velocity_i];
    vy_l   = prim_l[prim_lut.velocity_j];
    vz_l   = prim_l[prim_lut.velocity_k];
    Bx_l   = prim_l[prim_lut.bfield_i];
    By_l   = prim_l[prim_lut.bfield_j];
    Bz_l   = prim_l[prim_lut.bfield_k];

    Bv_l = Bx_l * vx_l + By_l * vy_l + Bz_l * vz_l;
    etot_l = cons_l[cons_lut.total_energy];
    pt_l = p_l + this->mag_pressure_(prim_l, prim_lut);
    cf_l = this->fast_magnetosonic_speed_(prim_l, prim_lut, eos);

    // load wr and compute the fast magnetosonic speed
    rho_r   = prim_r[prim_lut.density];
    p_r     = prim_r[prim_lut.pressure];
    vx_r    = prim_r[prim_lut.velocity_i];
    vy_r    = prim_r[prim_lut.velocity_j];
    vz_r    = prim_r[prim_lut.velocity_k];
    Bx_r    = prim_r[prim_lut.bfield_i];
    By_r    = prim_r[prim_lut.bfield_j];
    Bz_r    = prim_r[prim_lut.bfield_k];

    Bv_r = Bx_r * vx_r + By_r * vy_r + Bz_r * vz_r;
    etot_r = cons_r[cons_lut.total_energy];
    pt_r = p_r + this->mag_pressure_(prim_r, prim_lut);
    cf_r = this->fast_magnetosonic_speed_(prim_r, prim_lut, eos);

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

    if (S_l > 0) {
      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = flux_l[field];
      }
      return;
    } else if (S_r < 0) {
      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = flux_r[field];
      }
      return;
    } 

    // next, the middle (contact) wave
    S_M = ((S_r - vx_r)*rho_r*vx_r - (S_l - vx_l)*
	   rho_l*vx_l - pt_r + pt_l)/((S_r - vx_r)*
				      rho_r - (S_l - vx_l)*rho_l);

    // finally, the intermediate (Alfven) waves

    rho_ls = rho_l * (S_l - vx_l)/(S_l - S_M);
    rho_rs = rho_r * (S_r - vx_r)/(S_r - S_M);

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
      setup_cons_ast_(Us, cons_lut, S_M, rho_ls, vy_ls, vz_ls, etot_ls, Bx,
		      By_ls, Bz_ls);

      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = \
	  flux_l[field] + S_l*(Us[field] - cons_l[field]);
      }
      return;
    } else if (S_rs <= 0 && S_r >= 0) {
      // USE F_rs
      setup_cons_ast_(Us, cons_lut, S_M, rho_rs, vy_rs, vz_rs, etot_rs, Bx,
		      By_rs, Bz_rs);

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
      setup_cons_ast_(Us, cons_lut, S_M, rho_ls, vy_ls, vz_ls, etot_ls, Bx,
		      By_ls, Bz_ls);
      setup_cons_ast_(Uss, cons_lut, S_M, rho_ls, vy_ss, vz_ss, etot_lss, Bx,
		      By_ss, Bz_ss);

      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = \
	  (flux_l[field] + S_ls*Uss[field] - (S_ls - S_l)*Us[field] -
	   S_l*cons_l[field]);
      }
      return;
    } else if (S_M <= 0 && S_rs >= 0) {
      // USE F_rss
      setup_cons_ast_(Us, cons_lut, S_M, rho_rs, vy_rs, vz_rs, etot_rs, Bx,
		      By_rs, Bz_rs);
      setup_cons_ast_(Uss, cons_lut, S_M, rho_rs, vy_ss, vz_ss, etot_rss, Bx,
		      By_ss, Bz_ss);

      for (int field = 0; field<n_keys; field++){
	flux_arrays[field](iz,iy,ix) = \
	  (flux_r[field] + S_rs*Uss[field] - (S_rs - S_r)*Us[field] -
	   S_r*cons_r[field]);
      }
      return;
    }
  }
  
private:
  inline void setup_cons_ast_(enzo_float cons[], const field_lut cons_lut,
			      const enzo_float speed, const enzo_float rho,
			      const enzo_float vy, const enzo_float vz,
			      const enzo_float etot, const enzo_float Bx,
			      const enzo_float By, const enzo_float Bz)
  {

    // Helper function that factors out the filling of the of the asterisked
    // and double asterisked conserved quantities
    // This function intentially omits the following:
    // if (DualEnergyFormalism){
    //   The value we need to assign to cons["internal_energy"] is not clear
    //   In the original Enzo code, they set it equal to:
    //    - rho_ls * eint_ls
    //    - rho_rs * eint_rs
    //    - rho_ls * eint_lss
    //    - rho_rs * eint_rss
    //   Although eint_ls, eint_rs, eint_lss, and eint_rss are all declared as
    //   local variables, they are never actually initialized
    // }
    cons[cons_lut.density] = rho;
    cons[cons_lut.momentum_i] = rho * speed;
    cons[cons_lut.momentum_j] = rho * vy;
    cons[cons_lut.momentum_k] = rho * vz;
    cons[cons_lut.total_energy] = etot;
    cons[cons_lut.bfield_i] = Bx;
    cons[cons_lut.bfield_j] = By;
    cons[cons_lut.bfield_k] = Bz;
  }


private:
  // Below are two float maps used in the calculation
  enzo_float *Us;
  enzo_float *Uss;
};

#endif /* ENZO_ENZO_RIEMANN_HLLD_HPP */
