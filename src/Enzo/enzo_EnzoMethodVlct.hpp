
#ifndef ENZO_ENZO_METHOD_VLCT_HPP
#define ENZO_ENZO_METHOD_VLCT_HPP
class EnzoMethodVlct : public Method {

  /// @class    EnzoMethodVlct
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate VL + CT MHD method

public: // interface

  /// Create a new EnzoMethodVlct object
  EnzoMethodVlct();

  //EnzoMethodVlct()
  //  : Method(),
  //    eos_(NULL),
  //    half_dt_recon_(NULL),
  //    full_dt_recon_(NULL),
  //    riemann_solver_(NULL)
  //{ }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodVlct);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodVlct (CkMigrateMessage *m)
    : Method (m),
      eos_(NULL),
      half_dt_recon_(NULL),
      full_dt_recon_(NULL),
      riemann_solver_(NULL)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "vlct"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // methods

  // not sure if I will pass field_ids and blocks or arrays
  // not sure if this should be static
  void compute_flux_(Block *block, int dim, std::vector<int> &prim_ids,
		     std::vector<int> &bface_ids, std::vector<int> &priml_ids,
		     std::vector<int> &primr_ids, std::vector<int> &flux_ids,
		     EnzoReconstructor &reconstructor);

  // compute the Electric fields using the fluxes and cell-centered
  // primitives
  void compute_efields_(Block *block, std::vector<int> &xflux_ids,
			std::vector<int> &yflux_ids,
			std::vector<int> &zflux_ids,
			int center_efield_id,
			std::vector<int> &efield_ids,
			std::vector<int> &prim_ids,
			EnzoConstrainedTransport &ct);

  // update the tracked quantities
  // cur_cons_ids indicates the values of the conserved ids to add flux to
  // out_cons_ids indicates the place to update the conserved ids (this can be
  // identical to cur_cons_ids)
  void update_quantities_(Block *block, std::vector<int> &xflux_ids,
			  std::vector<int> &yflux_ids,
			  std::vector<int> &zflux_ids,
			  std::vector<int> &efield_ids,
			  std::vector<int> &cur_cons_ids,
			  std::vector<int> &out_cons_ids,
			  std::vector<int> &cur_bface_ids,
			  std::vector<int> &out_bface_ids,
			  double dt);

  // allocate the temporary fields needed for scratch space and store their ids
  // efield_ids refer to efields centered on the edges of cells
  void allocate_temp_fields_(Block *block, std::vector<int> &prim_ids,
			     std::vector<int> &priml_ids,
			     std::vector<int> &primr_ids,
			     std::vector<int> &xflux_ids,
			     std::vector<int> &yflux_ids,
			     std::vector<int> &zflux_ids,
			     std::vector<int> &efield_ids,
			     int &center_efield_id,
			     std::vector<int> &temp_cons_ids,
			     std::vector<int> &temp_bface_ids);

  // deallocate the temporary fields used for scratch space
  void deallocate_temp_fields_(Block *block, std::vector<int> &prim_ids,
			       std::vector<int> &priml_ids,
			       std::vector<int> &primr_ids,
			       std::vector<int> &xflux_ids,
			       std::vector<int> &yflux_ids,
			       std::vector<int> &zflux_ids,
			       std::vector<int> &efield_ids,
			       int center_efield_id,
			       std::vector<int> &temp_cons_ids,
			       std::vector<int> &temp_bface_ids);

protected: // attributes

  EnzoEquationOfState *eos_;
  EnzoReconstructor *half_dt_recon_;
  EnzoReconstructor *full_dt_recon_;
  EnzoRiemannSolver *riemann_solver_;
  
};

#endif /* ENZO_ENZO_METHOD_VLCT_HPP */
