// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodM1ClosureRT.hpp
/// @author   William Hicks (whicks@ucsd.edu) 
/// @date     Mon Aug  16 16:14:38 PDT 2021
/// @brief    [\ref Enzo] Declaration of EnzoMethodM1Closure
///           Radiative transfer using M1 closure method as implemented
///           in the RAMSES-RT code (https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R/abstract)

#ifndef ENZO_ENZO_METHOD_M1_CLOSURE
#define ENZO_ENZO_METHOD_M1_CLOSURE

class M1Tables {

// The M1Tables class contains methods for reading and accessing data tables
//   relevant to EnzoMethodM1Closure

public:
  M1Tables();

  int hll_table_f     (int i, int j) const throw() { return hll_table_f_[100*i+j]; }
  int hll_table_theta (int i, int j) const throw() { return hll_table_theta_[100*i+j]; }

  double hll_table_lambda_min (int i, int j) const throw() { return hll_table_lambda_min_[100*i+j]; }
  double hll_table_lambda_max (int i, int j) const throw() { return hll_table_lambda_max_[100*i+j]; }
  double hll_table_col3       (int i, int j) const throw() { return hll_table_col3_[100*i+j]; }
  double hll_table_col4       (int i, int j) const throw() { return hll_table_col4_[100*i+j]; }

private:
  void read_hll_eigenvalues(std::string hll_file) throw(); 
 
  std::vector<int> hll_table_f_, hll_table_theta_;
  std::vector<double> hll_table_lambda_min_, hll_table_lambda_max_;
  std::vector<double> hll_table_col3_, hll_table_col4_;
};

//-----------------------------------------------

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


 /// calls inject_photons(), sets the groups' mean cross sections & energies,
  /// then launches a refresh
  ///
  /// @note
  /// When the cross section calculation method is "vernier_average", this
  /// actually initiates a global reduction for computing the groups' mean
  /// cross sections & energies. The callback executed following the reduction
  /// subsequently initiates the callback.
  void call_inject_photons(EnzoBlock * enzo_block) throw();

  /// solves the transport equation & photoheating and photoionization rates
  ///
  /// This is invoked by the entry method right after the refresh-step and
  /// drives the rest of the calculation
  void call_solve_transport_eqn(EnzoBlock * enzo_block) throw();

  /// updates each cell in the integrated fields (e.g. "photon_density",
  /// "flux_x", "flux_y", "flux_z") with the sum of the values from each of the
  /// group fields.
  void sum_group_fields (EnzoBlock * enzo_block) throw();

  /// this is called after the reduction required by the "vernier_average"
  /// calculation method (related to computing the average photon energy,
  /// recombination cross-sections, and energy cross-sections for each photon
  /// group)
  ///
  /// the photon-injection refresh starts after this method completes. This
  /// method is not invoked by the other approach for computing cross sections
  /// & energies
  void set_global_averages (EnzoBlock * enzo_block, CkReductionMsg * msg) throw();

protected: // methods
 
  /// numerically integrate f from a to b with Simpson's rule
  ///
  /// @param a,b bounds of the integral
  /// @param n the number of subintervals to use (must be positive)
  /// @param f the function to be integrated
  /// @param v1,v2,v3 the trailing arguments that are passed to f
  double integrate_simpson(double a, double b, int n, 
               std::function<double(double,double,double,int)> f,double v1,double v2,int v3) throw();

  /// compute the planck function.
  ///
  /// @param nu frequency
  /// @param T temperature
  /// @param clight the speed of light
  /// @param dependent_variable denotes the numerator. When this is 0, the
  ///     numerator is 1. When this is 1, it's appropriate for photon density,
  ///     When this is 2, it's appropriate for energy density
  double planck_function(double nu, double T, double clight, int dependent_variable) throw();

  /// gives the name of the ScalarData that holds the "average" photoionization
  /// cross-section of photon group i, appropriate for gas species j
  const std::string sigN_string(int i, int j) throw()
    {return "sigN_" + std::to_string(i) + std::to_string(j);}
  
  /// gives the name of the ScalarData that holds the "average" energy
  /// cross-section of photon group i, appropriate for gas species j
  const std::string sigE_string(int i, int j) throw()
    {return "sigE_" + std::to_string(i) + std::to_string(j);}

  /// gives the name of the ScalarData that holds the "average" energy of
  /// photons in photon group i
  const std::string eps_string (int i) throw()
    {return "eps_"  + std::to_string(i);}

  /// gives the name of the ScalarData that total mass*luminosity emitted by
  /// star particles (the luminosity is for photon group i)
  const std::string mL_string  (int i) throw()
    {return "mL_" + std::to_string(i);}

  
  //--------- INJECTION STEP -------

  /// get stellar temperature for a given stellar mass
  double get_star_temperature(double M) throw();

  /// computes the number-density of photons (in a give group) emitted by a
  /// star-particle of mass `pmass` and luminosity `plum` based on the
  /// user-specified stellar SED
  double get_radiation_custom(EnzoBlock * enzo_block,
             double energy, double pmass, double plum, 
             double dt, double inv_vol, int igroup) throw();

  /// computes the number-density of photons (in a give group) that emitted by
  /// a star-particle of mass `pmass` (assuming a black-body spectrum)
  double get_radiation_blackbody(EnzoBlock * enzo_block, double pmass, 
             double freq_lower, double freq_upper, double clight, 
             double dt, double cell_volume, int igroup) throw();

  /// computes the number density of star particles and deposits them on the
  /// grid
  void inject_photons(EnzoBlock * enzo_block, int igroup) throw();

  //--------- TRANSPORT STEP --------


  double flux_function (double U_l, double U_lplus1,
                        double Q_l, double Q_lplus1,
                        double clight, double lmin, double lmax, std::string type) throw();

  void compute_hll_eigenvalues(double f, double theta, double * lmin, double * lmax, double clight) throw();

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

  /// computes the photon-loss term from attenuation by local gas
  double D_add_attenuation ( EnzoBlock * enzo_block, double clight, int i, int igroup) throw(); 

  /// helper function used in C_add_recombination
  double get_alpha (double T, int species, char rec_case) throw();

  /// helper function used in C_add_recombination
  int get_b_boolean (double E_lower, double E_upper, int species) throw();

  /// computes the photon-creation term from recombination in local gas
  double C_add_recombination (EnzoBlock * enzo_block, 
                              enzo_float * T, int i, int igroup, double E_lower, double E_uppe) throw();

  /// Computes the photoionization cross-section of particles in a given gas
  /// species (specified by type) and for photons of energy E
  double sigma_vernier (double energy, int type) throw();

  /// Calculates photoionization and heating rates in each cell
  void get_photoionization_and_heating_rates (EnzoBlock * enzo_block, double clight) throw();

protected: // attributes
  int N_groups_;

  // Refresh id's
  int ir_injection_;

  // Tables relevant to M1 closure method
  M1Tables * M1_tables;
};


#endif /* ENZO_ENZO_METHOD_M1_CLOSURE_HPP */
