// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodRadiationInjection.hpp
/// @author   William Hicks (whicks@ucsd.edu) 
/// @date     Mon Aug  16 16:14:38 PDT 2021
/// @brief    [\ref Enzo] Declaration of EnzoMethodRadiationInjection
///           Radiative transfer using M1 closure method as implemented
///           in the RAMSES-RT code

#ifndef ENZO_ENZO_METHOD_RADIATION_INJECTION
#define ENZO_ENZO_METHOD_RADIATION_INJECTION

class EnzoMethodRadiationInjection : public Method {

  /// @class    EnzoMethodRadiationInjection 
  /// @ingroup  Enzo
  ///
  /// @brief    [\ref Enzo] Declaration of EnzoMethodRadiationInjection
  ///           Radiative transfer using M1 closure method as implemented
  ///           in the RAMSES-RT code. This is the photon injection step

public: // interface

  /// Create a new EnzoMethodRadiationInjection object

  EnzoMethodRadiationInjection();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodRadiationInjection );
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodRadiationInjection  (CkMigrateMessage *m)
    : Method (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "radiation_injection"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block ) const throw();

protected: // methods
  void compute_ (Block * block) throw();

protected: // attributes

};

#endif /* ENZO_ENZO_METHOD_RADIATION_INJECTION_HPP */
