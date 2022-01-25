// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRamsesRT.hpp
/// @author   William Hicks (whicks@ucsd.edu) 
/// @date     Mon Aug  16 16:14:38 PDT 2021
/// @brief    [\ref Enzo] Declaration of EnzoMethodRamsesRT
///           Radiative transfer using M1 closure method as implemented
///           in the RAMSES-RT code

#ifndef ENZO_ENZO_METHOD_RAMSES_RT
#define ENZO_ENZO_METHOD_RAMSES_RT

class EnzoMethodRamsesRT : public Method {

  /// @class    EnzoMethodRamsesRT 
  /// @ingroup  Enzo
  ///
  /// @brief    [\ref Enzo] Declaration of EnzoMethodRamsesRT
  ///           Radiative transfer using M1 closure method as implemented
  ///           in the RAMSES-RT code

public: // interface

  /// Create a new EnzoMethodRamsesRT object

  EnzoMethodRamsesRT(const int N_groups, const double clight);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodRamsesRT );
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodRamsesRT  (CkMigrateMessage *m)
    : Method (m)
    , N_groups_(0)
    , clight_(0.0)
    , eps_()
    , sigN_() 
    , sigE_()
    , gfracN_()
    , gfracF_()
    , ir_injection_(-1)
    , ir_transport_(-1)
   // , igroup_()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "ramses_rt"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block ) const throw();


  //--------- CONTROL FLOW --------

  // calls inject_photons(), then does a refresh
  void call_inject_photons(EnzoBlock * enzo_block) throw();

  void call_solve_transport_eqn(EnzoBlock * enzo_block) throw();

  void sum_group_fields (EnzoBlock * enzo_block) throw();

protected: // methods

  double integrate_simpson(double a, double b, int n, 
               std::function<double(double,double,double,int)> f,double v1,double v2,int v3) throw();

  double planck_function(double nu, double T, double clight, int dependent_variable) throw();

  //--------- INJECTION STEP -------

  double get_star_temperature(double M) throw();

  void get_radiation_blackbody(EnzoBlock * enzo_block, enzo_float * N, int i, double T, 
             double freq_lower, double freq_upper, double clight, double f_esc) throw();

  void inject_photons(EnzoBlock * enzo_block) throw();


  //--------- TRANSPORT STEP --------


  double flux_function (double U_l, double U_lplus1,
                        double Q_l, double Q_lplus1,
                        double clight, std::string type) throw();

  double deltaQ_faces (double U_l, double U_lplus1, double U_lminus1, 
                       double Q_l, double Q_lplus1, double Q_lminus1,
                       double clight) throw();

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

  void solve_transport_eqn (EnzoBlock * enzo_block) throw();

  //---------- THERMOCHEMISTRY STEP ------------
  // Interaction with matter is completely local, so don't need a refresh before this step
  void add_attenuation ( EnzoBlock * enzo_block, 
              enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz, 
              double clight, int i) throw(); 

  double get_beta (double T, int species) throw();

  double get_alpha (double T, int species, char rec_case) throw();

  int get_b_boolean (double E_lower, double E_upper, int species) throw();

  void recombination_photons (EnzoBlock * enzo_block, enzo_float * N, 
                              enzo_float * T, int i, double E_lower, double E_upper) throw();

  void recombination_chemistry (EnzoBlock * enzo_block) throw();

  double sigma_vernier (double energy, int type) throw();

  void get_photoionization_and_heating_rates (EnzoBlock * enzo_block, double clight) throw();

  void compute_ (Block * block) throw();


protected: // attributes
  int N_groups_;
  double clight_;
  std::vector<double> eps_, sigN_, sigE_;
  std::vector<double> gfracN_, gfracF_;

  // Refresh id's
  int ir_injection_, ir_transport_;
};

#endif /* ENZO_ENZO_METHOD_RAMSES_RT_HPP */
