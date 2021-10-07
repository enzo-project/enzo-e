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

protected: // methods

  //--------- INJECTION STEP -------

  double get_star_temperature(double M) throw();

  void get_radiation_blackbody(enzo_float * N, int i, double T, 
             double freq_lower, double freq_upper, double dt, double clight, double f_esc) throw();

  void inject_photons(EnzoBlock * enzo_block) throw();


  //--------- TRANSPORT STEP --------


  double flux_function (double U_l, double U_lplus1,
                        double Q_l, double Q_lplus1,
                        std::string type) throw();

  double deltaQ_faces (double U_l, double U_lplus1, double U_lminus1, 
                       double Q_l, double Q_lplus1, double Q_lminus1)
                       throw();

  void get_reduced_variables (double * chi_idx, double (*n_idx)[3], int i, double clight,
                              enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz) 
                              throw(); 

  void get_pressure_tensor (EnzoBlock * enzo_block, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz)
                       throw();

  void get_U_update (EnzoBlock * enzo_block, double * N_update, 
                       double * Fx_update, double * Fy_update, double * Fz_update, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz,
                       double hx, double hy, double hz, double dt, double clight, 
                       int i, int idx, int idy, int idz) throw();

  void solve_transport_eqn (EnzoBlock * enzo_block) throw();

  void compute_ (Block * block) throw();
protected: // attributes
  int N_groups_;
  //int igroup_;// = 0;
  double clight_;
  std::vector<double> eps_, sigN_, sigE_;
  std::vector<double> gfracN_, gfracF_;

  // Refresh id's
  int ir_injection_, ir_transport_;
};

#endif /* ENZO_ENZO_METHOD_RAMSES_RT_HPP */
