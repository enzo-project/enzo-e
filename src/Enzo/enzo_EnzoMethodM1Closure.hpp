// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodM1ClosureRT.hpp
/// @author   William Hicks (whicks@ucsd.edu) 
/// @date     Mon Aug  16 16:14:38 PDT 2021
/// @brief    [\ref Enzo] Declaration of EnzoMethodM1Closure
///           Radiative transfer using M1 closure method as implemented
///           in the RAMSES-RT code (https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R/abstract)

#ifndef ENZO_ENZO_METHOD_M1_CLOSURE
#define ENZO_ENZO_METHOD_M1_CLOSURE

class EnzoMethodM1Closure : public Method {

  /// @class    EnzoMethodM1Closure
  /// @ingroup  Enzo
  ///
  /// @brief    [\ref Enzo] Declaration of EnzoMethodM1Closure
  ///           Radiative transfer using M1 closure method as implemented
  ///           in the RAMSES-RT code

public: // interface

  /// Create a new EnzoMethodM1Closure object

  EnzoMethodM1Closure(const int N_groups);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodM1Closure );
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodM1Closure (CkMigrateMessage *m)
    : Method (m)
    , N_groups_(0)
    , ir_injection_(-1)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "m1_closure"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block ) throw();


  //--------- CONTROL FLOW --------
  //  compute_ -> call_inject_photons -> inject_photons ->
  //  refresh -> call_solve_transport_eqn -> 
  //  solve_transport_eqn, C_add_recombination, D_add_attenuation, 
  //    get_photoionization_and_heating_rates 


  // calls inject_photons(), then does a refresh
  void call_inject_photons(EnzoBlock * enzo_block) throw();

  void call_solve_transport_eqn(EnzoBlock * enzo_block) throw();

  void sum_group_fields (EnzoBlock * enzo_block) throw();

  void set_global_averages (EnzoBlock * enzo_block, CkReductionMsg * msg) throw();

protected: // methods
 
  double integrate_simpson(double a, double b, int n, 
               std::function<double(double,double,double,int)> f,double v1,double v2,int v3) throw();

  double planck_function(double nu, double T, double clight, int dependent_variable) throw();

  const std::string sigN_string(int i, int j) throw()
    {return "sigN_" + std::to_string(i) + std::to_string(j);}
  const std::string sigE_string(int i, int j) throw()
    {return "sigE_" + std::to_string(i) + std::to_string(j);}
  const std::string eps_string (int i) throw()
    {return "eps_"  + std::to_string(i);}
  const std::string mL_string  (int i) throw()
    {return "mL_" + std::to_string(i);}

  
  //--------- INJECTION STEP -------

  double get_star_temperature(double M) throw();

  double get_radiation_custom(EnzoBlock * enzo_block,
             double energy, double pmass, double plum, 
             double dt, double inv_vol, int i, int igroup) throw();

  double get_radiation_blackbody(EnzoBlock * enzo_block, double pmass, 
             double freq_lower, double freq_upper, double clight, 
             double dt, double cell_volume, int i, int igroup) throw();

  void inject_photons(EnzoBlock * enzo_block, int igroup) throw();

  //--------- TRANSPORT STEP --------


  double flux_function (double U_l, double U_lplus1,
                        double Q_l, double Q_lplus1,
                        double clight, double lmin, double lmax, std::string type) throw();

  double compute_hll_eigenvalues(double f, double theta, double * lmin, double * lmax, double clight) throw();

  double deltaQ_faces (double U_l, double U_lplus1, double U_lminus1, 
                       double Q_l, double Q_lplus1, double Q_lminus1,
                       double clight, double lmin, double lmax, std::string flux_type) throw();

  void get_reduced_variables (double * chi_idx, double (*n_idx)[3], int i, double clight,
                              enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz) 
                              throw(); 

  void get_pressure_tensor (EnzoBlock * enzo_block, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz,
                       double clight) throw();

  void get_U_update (EnzoBlock * enzo_block, double * N_update, 
                       double * Fx_update, double * Fy_update, double * Fz_update, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz,
                       double hx, double hy, double hz, double dt, double clight, 
                       int i, int idx, int idy, int idz) throw();

  void solve_transport_eqn (EnzoBlock * enzo_block, int igroup) throw();

  void add_LWB (EnzoBlock * enzo_block, double J21);

  //---------- THERMOCHEMISTRY STEP ------------


  // Interaction with matter is completely local, so don't need a refresh before this step
  void D_add_attenuation ( EnzoBlock * enzo_block, double * D, double clight, int i, int igroup) throw(); 

  double get_alpha (double T, int species, char rec_case) throw();

  int get_b_boolean (double E_lower, double E_upper, int species) throw();

  void C_add_recombination (EnzoBlock * enzo_block, double * C, 
                              enzo_float * T, int i, int igroup, double E_lower, double E_uppe) throw();

  double sigma_vernier (double energy, int type) throw();

  void get_photoionization_and_heating_rates (EnzoBlock * enzo_block, double clight) throw();

  //--------------------------

  void compute_ (Block * block) throw();


protected: // attributes
  int N_groups_;

  // Refresh id's
  int ir_injection_;
};

#endif /* ENZO_ENZO_METHOD_M1_CLOSURE_HPP */
