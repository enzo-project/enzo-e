#ifndef ENZO_ENZO_RIEMANN_HLLD_HPP
#define ENZO_ENZO_RIEMANN_HLLD_HPP

// This is a quick port from Enzo's Riemann_HLLD_MHD.C, written by J. S. Oishi
// It needs to be rewritten to avoid duplicate calculations

class EnzoRiemannHLLD : public EnzoRiemann
{

  /// @class    EnzoRiemannHLLD
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates HLLE approximate Riemann Solver

public: // interface
  /// Create a new EnzoRiemannHLLE object
  EnzoRiemannHLLD() throw()
  : EnzoRiemann()
  { }

  EnzoRiemannHLLD(std::vector<std::string> &extra_scalar_groups,
		  std::vector<std::string> &extra_vector_groups,
		  std::vector<std::string> &extra_passive_groups,
		  FluxFunctor** flux_funcs, int n_funcs)
    : EnzoRiemann(extra_scalar_groups, extra_vector_groups,
		  extra_passive_groups, flux_funcs, n_funcs)
  { }

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
  };

  void calc_riemann_fluxes_(flt_map &flux_l, flt_map &flux_r,
			    flt_map &prim_l, flt_map &prim_r,
			    flt_map &cons_l, flt_map &cons_r,
			    std::vector<std::string> &cons_keys,
			    std::size_t n_keys,
			    EnzoEquationOfState *eos,
			    const int iz, const int iy, const int ix,
			    array_map &flux_arrays)
  {
    // Need to edit so that we can make the flt_maps const

    // NEED TO DEFINE BFLOAT_EPSILON
    // I think this determines how close a value needs to be to 0 to count as 0
    // In that case, it may be worthwhile to use it to determine the upwind
    // direction
    enzo_float BFLOAT_EPSILON = 0;
    // we should either store Us and Uss as members of the data structure OR, we
    // should pass them in
    flt_map Us, Uss;
    Us.reserve(n_keys); Uss.reserve(n_keys);

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
    bool DualEnergyFormalism = false;

    // First, compute Fl and Ul
    rho_l  = prim_l["density"];
    if (DualEnergyFormalism){
      //eint_l = prim_l["internal_energy"];
      //enzo_float h, dpdrho, dpde;
      //EOS(p_l, rho_l, eint_l, h, cs_l, dpdrho, dpde, EOSType, 2);
    } else {
      p_l = prim_l["pressure"];
      //cs_l = this->sound_speed_(prim_l, eos);
    }
    vx_l   = prim_l["velocity_i"];
    vy_l   = prim_l["velocity_j"];
    vz_l   = prim_l["velocity_k"];    
    Bx_l   = prim_l["bfield_i"];
    By_l   = prim_l["bfield_j"];
    Bz_l   = prim_l["bfield_k"];

    Bv_l = Bx_l * vx_l + By_l * vy_l + Bz_l * vz_l;
    etot_l = cons_l["total_energy"];
    pt_l = p_l + this->mag_pressure_(prim_l);
    cf_l = this->fast_magnetosonic_speed_(prim_l, eos);

    //B2 = Bx_l * Bx_l + By_l * By_l + Bz_l * Bz_l;
    //v2 = vx_l * vx_l + vy_l * vy_l + vz_l * vz_l;
    //etot_l = rho_l * (eint_l + 0.5 * v2) + 0.5 * B2;
    //cf_l = std::sqrt((gamma * p_l + B2 +
    //		 std::sqrt((gamma * p_l + B2) * (gamma * p_l + B2)
    //		      - 4. * gamma * p_l * Bx_l * Bx_l))/(2. * rho_l));
    

    //cons_l["density"] = rho_l;
    //cons_l["momentum_i"] = rho_l * vx_l;
    //cons_l["momentum_j"] = rho_l * vy_l;
    //cons_l["momentum_k"] = rho_l * vz_l;
    //cons_l["total_energy"] = etot_l;
    //if (DualEnergyFormalism) {
    //  cons_l["internal_energy"] = rho_l * eint_l;
    //}
    //cons_l["bfield_i"] = Bx_l;
    //cons_l["bfield_j"] = By_l;
    //cons_l["bfield_k"] = Bz_l;

    // The following calculations are redundant (they should have been handled
    // prior to the function call - this is for debugging)
    flux_l["density"] = rho_l * vx_l;
    flux_l["momentum_i"] = cons_l["momentum_i"] * vx_l + pt_l - Bx_l * Bx_l;
    flux_l["momentum_j"] = cons_l["momentum_j"] * vx_l - Bx_l * By_l;
    flux_l["momentum_k"] = cons_l["momentum_k"] * vx_l - Bx_l * Bz_l;
    flux_l["total_energy"] = (etot_l + pt_l) * vx_l - Bx_l * Bv_l;
    if (DualEnergyFormalism) {
      flux_l["internal_energy"] = cons_l["internal_energy"] * vx_l;
    }
    flux_l["bfield_i"] = 0.0;
    flux_l["bfield_j"] = vx_l*By_l - vy_l*Bx_l;
    flux_l["bfield_k"] = -vz_l*Bx_l + vx_l*Bz_l;

    //compute Ur and Fr
    rho_r   = prim_r["density"];
    if (DualEnergyFormalism){
      //eint_r = prim_r["internal_energy"];
      //enzo_float h, dpdrho, dpde;
      //EOS(p_l, rho_l, eint_l, h, cs_l, dpdrho, dpde, EOSType, 2);
    } else {
      p_r = prim_r["pressure"];
      //cs_r = this->sound_speed_(prim_r, eos);
    }
    vx_r    = prim_r["velocity_i"];
    vy_r    = prim_r["velocity_j"];
    vz_r    = prim_r["velocity_k"];
    Bx_r    = prim_r["bfield_i"];
    By_r    = prim_r["bfield_j"];
    Bz_r    = prim_r["bfield_k"];

    Bv_r = Bx_r * vx_r + By_r * vy_r + Bz_r * vz_r;
    etot_r = cons_r["total_energy"];
    pt_r = p_r + this->mag_pressure_(prim_r);
    cf_r = this->fast_magnetosonic_speed_(prim_r,eos);

    // The following calculations are redundant (they should have been handled
    // prior to the function call - this is for debugging)
    flux_r["density"] = rho_r * vx_r;
    flux_r["momentum_i"] = cons_r["momentum_i"] * vx_r + pt_r - Bx_r * Bx_r;
    flux_r["momentum_j"] = cons_r["momentum_j"] * vx_r - Bx_r * By_r;
    flux_r["momentum_k"] = cons_r["momentum_k"] * vx_r - Bx_r * Bz_r;
    flux_r["total_energy"] = (etot_r + pt_r) * vx_r - Bx_r * Bv_r;
    if (DualEnergyFormalism) {
      flux_r["internal_energy"] = cons_r["internal_energy"] * vx_l;
    }
    flux_r["bfield_i"] = 0.0;
    flux_r["bfield_j"] = vx_r * By_r - vy_r * Bx_r;
    flux_r["bfield_k"] = -vz_r * Bx_r + vx_r * Bz_r;

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
	flux_arrays[key](iz,iy,ix) = flux_l[key];
      }
      return;
    } else if (S_r < 0) {
      for (std::size_t field_ind = 0; field_ind<n_keys; field_ind++){
	std::string key = cons_keys[field_ind];
	flux_arrays[key](iz,iy,ix) = flux_r[key];
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
	  flux_l[key] + S_l*(Us[key] - cons_l[key]);
      }
      return;
    } else if (S_rs <= 0 && S_r >= 0) {
      // USE F_rs
      setup_cons_ast_(Us, S_M, rho_rs, vy_rs, vz_rs, etot_rs, Bx, By_rs, Bz_rs);

      for (std::size_t field_ind = 0; field_ind<n_keys; field_ind++){
	std::string key = cons_keys[field_ind];
	flux_arrays[key](iz,iy,ix) = \
	  flux_r[key] + S_r*(Us[key] - cons_r[key]);
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
	  flux_l[key] + S_ls*Uss[key] - (S_ls - S_l)*Us[key] - S_l*cons_l[key];
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
	  flux_r[key] + S_rs*Uss[key] - (S_rs - S_r)*Us[key] - S_r*cons_r[key];
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
		       
};

#endif /* ENZO_ENZO_RIEMANN_HLLD_HPP */
