// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPpm.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo PPM hydro method

#ifndef ENZO_ENZO_METHOD_PPM_HPP
#define ENZO_ENZO_METHOD_PPM_HPP

class EnzoMethodPpm : public Method {

  /// @class    EnzoMethodPpm
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's PPM hydro method

public: // interface

  /// Create a new EnzoMethodPpm object
  EnzoMethodPpm(bool store_fluxes_for_corrections,
                bool diffusion,
                int flattening,
                bool pressure_free,
                bool steepening,
                bool use_minimum_pressure_support,
                enzo_float minimum_pressure_support_parameter);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodPpm);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodPpm (CkMigrateMessage *m)
    : Method (m),
      comoving_coordinates_(false),
      store_fluxes_for_corrections_(false),
      diffusion_(false),
      flattening_(0),
      pressure_free_(false),
      steepening_(false),
      use_minimum_pressure_support_(false),
      minimum_pressure_support_parameter_(0.0)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual methods

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "ppm"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();

protected: // interface

  bool comoving_coordinates_;
  bool store_fluxes_for_corrections_;

  bool diffusion_;
  int flattening_;
  bool pressure_free_;
  bool steepening_;

  // in the future, it may make sense to move the following parameters to
  // PhysicsFluidProps or EnzoFluidFloorsConfig
  bool use_minimum_pressure_support_;
  enzo_float minimum_pressure_support_parameter_;
};

#endif /* ENZO_ENZO_METHOD_PPM_HPP */
