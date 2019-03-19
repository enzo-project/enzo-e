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
  /// Create a new EnzoRiemannHLLE object
  EnzoRiemannHLLD() throw()
  : EnzoRiemann()
  {
    // Reserve space in the following scratch space flt_maps
    Us.reserve(n_keys_);
    Uss.reserve(n_keys_);
  }

  EnzoRiemannHLLD(std::vector<std::string> &extra_scalar_groups,
		  std::vector<std::string> &extra_vector_groups,
		  std::vector<std::string> &extra_passive_groups,
		  FluxFunctor** flux_funcs, int n_funcs)
    : EnzoRiemann(extra_scalar_groups, extra_vector_groups,
		  extra_passive_groups, flux_funcs, n_funcs)
  {
    // Reserve space in the following scratch space flt_maps
    Us.reserve(n_keys_);
    Uss.reserve(n_keys_);
  }

  /// Virtual destructor
  virtual ~EnzoRiemannHLLD()
  {  }

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
    // Reserve space in the following scratch space flt_maps
    Us.reserve(n_keys_);
    Uss.reserve(n_keys_);
  };

  void calc_riemann_fluxes_(const flt_map &flux_l, const flt_map &flux_r,
			    const flt_map &prim_l, const flt_map &prim_r,
			    const flt_map &cons_l, const flt_map &cons_r,
			    std::vector<std::string> &cons_keys,
			    std::size_t n_keys,
			    EnzoEquationOfState *eos,
			    const int iz, const int iy, const int ix,
			    array_map &flux_arrays)
  {

    // NEED TO DEFINE BFLOAT_EPSILON
    // I think this determines how close a value needs to be to 0 to count as 0
    // In that case, it may be worthwhile to use it to determine the upwind
    // direction
    enzo_float BFLOAT_EPSILON = 0;
    // This method makes use of the member variables Us and Uss

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
    rho_l  = prim_l.at("density");
    p_l    = prim_l.at("pressure");
    vx_l   = prim_l.at("velocity_i");
    vy_l   = prim_l.at("velocity_j");
    vz_l   = prim_l.at("velocity_k");
    Bx_l   = prim_l.at("bfield_i");
    By_l   = prim_l.at("bfield_j");
    Bz_l   = prim_l.at("bfield_k");

    Bv_l = Bx_l * vx_l + By_l * vy_l + Bz_l * vz_l;
    etot_l = cons_l.at("total_energy");
    pt_l = p_l + this->mag_pressure_(prim_l);
    cf_l = this->fast_magnetosonic_speed_(prim_l, eos);

    // load wr and compute the fast magnetosonic speed
    rho_r   = prim_r.at("density");
    p_r     = prim_r.at("pressure");
    vx_r    = prim_r.at("velocity_i");
    vy_r    = prim_r.at("velocity_j");
    vz_r    = prim_r.at("velocity_k");
    Bx_r    = prim_r.at("bfield_i");
    By_r    = prim_r.at("bfield_j");
    Bz_r    = prim_r.at("bfield_k");

    Bv_r = Bx_r * vx_r + By_r * vy_r + Bz_r * vz_r;
    etot_r = cons_r.at("total_energy");
    pt_r = p_r + this->mag_pressure_(prim_r);
    cf_r = this->fast_magnetosonic_speed_(prim_r,eos);

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
      for (std::size_t field_ind = 0; field_ind<n_keys; field_ind++){
	std::string key = cons_keys[field_ind];
	flux_arrays[key](iz,iy,ix) = flux_l.at(key);
      }
      return;
    } else if (S_r < 0) {
      for (std::size_t field_ind = 0; field_ind<n_keys; field_ind++){
	std::string key = cons_keys[field_ind];
	flux_arrays[key](iz,iy,ix) = flux_r.at(key);
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

    S_ls = S_M - std::abs(Bx)/std::sqrt(rho_ls);
    S_rs = S_M + std::abs(Bx)/std::sqrt(rho_rs);

    pt_s = ((S_r -  vx_r) * rho_r*pt_l - (S_l - vx_l) * rho_l * pt_r +
	    rho_l*rho_r*(S_r - vx_r)*(S_l - vx_l)*
	    (vx_r - vx_l))/((S_r - vx_r)*rho_r - (S_l - vx_l)*rho_l);

    sam = vx_l - cf_l;
    sap = vx_l + cf_l;
      
    if ((std::abs(S_M - vx_l) <= BFLOAT_EPSILON) and 
        (std::abs(By_l) <= BFLOAT_EPSILON) and 
        (std::abs(Bz_l) <= BFLOAT_EPSILON) and 
        (Bx*Bx >= gamma * p_l) and
        ((std::abs(S_l - sam) <= BFLOAT_EPSILON) or
	 (std::abs(S_l - sap) <= BFLOAT_EPSILON)) ) {
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
      
    if ((std::abs(S_M - vx_r) <= BFLOAT_EPSILON) and 
        (std::abs(By_r) <= BFLOAT_EPSILON) and 
        (std::abs(Bz_r) <= BFLOAT_EPSILON) and 
        (Bx*Bx >= gamma * p_r) and
        ((std::abs(S_r - sam) <= BFLOAT_EPSILON) or
	 (std::abs(S_r - sap) <= BFLOAT_EPSILON)) ) {
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
      setup_cons_ast_(Us, S_M, rho_ls, vy_ls, vz_ls, etot_ls, Bx, By_ls, Bz_ls);

      for (std::size_t field_ind = 0; field_ind<n_keys; field_ind++){
	std::string key = cons_keys[field_ind];
	flux_arrays[key](iz,iy,ix) = \
	  flux_l.at(key) + S_l*(Us[key] - cons_l.at(key));
      }
      return;
    } else if (S_rs <= 0 && S_r >= 0) {
      // USE F_rs
      setup_cons_ast_(Us, S_M, rho_rs, vy_rs, vz_rs, etot_rs, Bx, By_rs, Bz_rs);

      for (std::size_t field_ind = 0; field_ind<n_keys; field_ind++){
	std::string key = cons_keys[field_ind];
	flux_arrays[key](iz,iy,ix) = \
	  flux_r.at(key) + S_r*(Us[key] - cons_r.at(key));
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
      setup_cons_ast_(Us, S_M, rho_ls, vy_ls, vz_ls, etot_ls, Bx, By_ls, Bz_ls);
      setup_cons_ast_(Uss, S_M, rho_ls, vy_ss, vz_ss, etot_lss, Bx,
		      By_ss, Bz_ss);

      for (std::size_t field_ind = 0; field_ind<n_keys; field_ind++){
	std::string key = cons_keys[field_ind];
	flux_arrays[key](iz,iy,ix) = \
	  (flux_l.at(key) + S_ls*Uss[key] - (S_ls - S_l)*Us[key] -
	   S_l*cons_l.at(key));
      }
      return;
    } else if (S_M <= 0 && S_rs >= 0) {
      // USE F_rss
      setup_cons_ast_(Us, S_M, rho_rs, vy_rs, vz_rs, etot_rs, Bx, By_rs, Bz_rs);
      setup_cons_ast_(Uss, S_M, rho_rs, vy_ss, vz_ss, etot_rss, Bx,
		      By_ss, Bz_ss);

      for (std::size_t field_ind = 0; field_ind<n_keys; field_ind++){
	std::string key = cons_keys[field_ind];
	flux_arrays[key](iz,iy,ix) = \
	  (flux_r.at(key) + S_rs*Uss[key] - (S_rs - S_r)*Us[key] -
	   S_r*cons_r.at(key));
      }
      return;
    }
  }
  
private:
  inline void setup_cons_ast_(flt_map& cons, const enzo_float speed,
			      const enzo_float rho, const enzo_float vy,
			      const enzo_float vz, const enzo_float etot,
			      const enzo_float Bx, const enzo_float By,
			      const enzo_float Bz)
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
    cons["density"] = rho;
    cons["momentum_i"] = rho * speed;
    cons["momentum_j"] = rho * vy;
    cons["momentum_k"] = rho * vz;
    cons["total_energy"] = etot;
    cons["bfield_i"] = Bx;
    cons["bfield_j"] = By;
    cons["bfield_k"] = Bz;
  }


private:
  // Below are two float maps used in the calculation
  flt_map Us;
  flt_map Uss;
};

#endif /* ENZO_ENZO_RIEMANN_HLLD_HPP */
