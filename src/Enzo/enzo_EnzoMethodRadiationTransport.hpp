// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRadiationTransport.hpp
/// @author   William Hicks (whicks@ucsd.edu) 
/// @date     Mon Aug  16 16:14:38 PDT 2021
/// @brief    [\ref Enzo] Declaration of EnzoMethodRadiationTransport
///           Radiative transfer using M1 closure method as implemented
///           in the RAMSES-RT code

#ifndef ENZO_ENZO_METHOD_RADIATION_TRANSPORT
#define ENZO_ENZO_METHOD_RADIATION_TRANSPORT

class EnzoMethodRadiationTransport : public Method {

  /// @class    EnzoMethodRadiationTransport 
  /// @ingroup  Enzo
  ///
  /// @brief    [\ref Enzo] Declaration of EnzoMethodRadiationTransport
  ///           Radiative transfer using M1 closure method as implemented
  ///           in the RAMSES-RT code

public: // interface

  /// Create a new EnzoMethodRadiationTransport object

  EnzoMethodRadiationTransport();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodRadiationTransport );
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodRadiationTransport  (CkMigrateMessage *m)
    : Method (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "radiation_transport"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block ) const throw();

protected: // methods


  //--------- TRANSPORT STEP --------


  double flux_function (double U_l, double U_lplus1,
                        double Q_l, double Q_lplus1,
                        std::string type, double clight) throw();

  double deltaQ_faces (double U_l, double U_lplus1, double U_lminus1, 
                       double Q_l, double Q_lplus1, double Q_lminus1,
                       double clight) throw();

  void get_reduced_variables (double * chi_idx, double (*n_idx)[3], int i, double clight,
                              enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz) 
                              throw(); 

  void get_pressure_tensor (Block * block, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz,
                       double clight) throw();

  void get_U_update (Block * block, double * N_update, 
                       double * Fx_update, double * Fy_update, double * Fz_update, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz, 
                       double hx, double hy, double hz, double dt, double clight, 
                       int i, int idx, int idy, int idz) throw();

  void transport_photons (Block * block, double clight) throw();

  void compute_ (Block * block) throw();
protected: // attributes

};

#endif /* ENZO_ENZO_METHOD_RADIATION_TRANSPORT_HPP */
