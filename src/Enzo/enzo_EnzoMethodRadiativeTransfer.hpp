// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRadiativeTransfer.hpp
/// @author   William Hicks (whicks@ucsd.edu) 
/// @date     Mon Aug  16 16:14:38 PDT 2021
/// @brief    [\ref Enzo] Declaration of EnzoMethodRadiativeTransfer
///           Radiative transfer using M1 closure method as implemented
///           in the RAMSES-RT code

#ifndef ENZO_ENZO_METHOD_RADIATIVE_TRANSFER
#define ENZO_ENZO_METHOD_RADIATIVE_TRANSFER

class EnzoMethodRadiativeTransfer : public Method {

  /// @class    EnzoMethodRadiativeTransfer 
  /// @ingroup  Enzo
  ///
  /// @brief    [\ref Enzo] Declaration of EnzoMethodRadiativeTransfer
  ///           Radiative transfer using M1 closure method as implemented
  ///           in the RAMSES-RT code

public: // interface

  /// Create a new EnzoMethodRadiativeTransfer object

  EnzoMethodRadiativeTransfer();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodRadiativeTransfer );
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodRadiativeTransfer  (CkMigrateMessage *m)
    : Method (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "radiative_transfer"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block ) const throw();

protected: // methods


  // -------- INJECTION STEP -----------
  void inject_photons (Block * block) throw();

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

  void get_updated_values (Block * block, double * N_update, 
                       double * Fx_update, double * Fy_update, double * Fz_update, 
                       enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz, 
                       double hx, double hy, double hz, double dt, double clight, 
                       int i, int idx, int idy, int idz) throw();


  void update_fluxes_1D (double * N_update, double * F_update, 
                         enzo_float * N, enzo_float * Fx, enzo_float * Fy, enzo_float * Fz,
                         double h, double dt, double clight, int i, int increment, int axis) throw(); 
  
  void transport_photons (Block * block, double clight) throw();
  // --------- CHEMISTRY STEP ---------

  void thermochemistry (Block * block) throw(); 



  void compute_ (Block * block) throw();
protected: // attributes

};

#endif /* ENZO_ENZO_METHOD_RADIATIVE_TRANSFER_HPP */
