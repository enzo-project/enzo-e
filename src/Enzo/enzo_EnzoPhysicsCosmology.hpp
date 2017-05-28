// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoPhysicsCosmology.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-05-23
/// @brief    [\ref Enzo] Declaration of the EnzoPhysicsCosmology class

#ifndef ENZO_ENZO_PHYSICS_COSMOLOGY_HPP
#define ENZO_ENZO_PHYSICS_COSMOLOGY_HPP

class EnzoPhysicsCosmology : public Physics {

  /// @class    EnzoPhysicsCosmology
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoPhysicsCosmology() throw()
  : Physics(),
    hubble_constant_now_(0.701),
    omega_matter_now_(0.279),
    omega_dark_matter_now_(-1.0),
    omega_lambda_now_(0.721),
    comoving_box_size_(64),
    max_expansion_rate_(0.01),
    initial_redshift_(20.0),
    final_redshift_(0.0)
  {  }

  /// Constructor
  EnzoPhysicsCosmology
  (
   enzo_float hubble_constant_now,
   enzo_float omega_matter_now,
   enzo_float omega_dark_matter_now,
   enzo_float omega_lambda_now,
   enzo_float comoving_box_size,
   enzo_float max_expansion_rate,
   enzo_float initial_redshift,
   enzo_float final_redshift
   )
    : Physics(),
      hubble_constant_now_(hubble_constant_now),
      omega_matter_now_(omega_matter_now),
      omega_dark_matter_now_(omega_dark_matter_now),
      omega_lambda_now_(omega_lambda_now),
      comoving_box_size_(comoving_box_size),
      max_expansion_rate_(max_expansion_rate),
      initial_redshift_(initial_redshift),
      final_redshift_(final_redshift)
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoPhysicsCosmology);

  /// CHARM++ migration constructor
  EnzoPhysicsCosmology(CkMigrateMessage *m)
    : Physics (m)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    Physics::pup(p);
    
    p | hubble_constant_now_;
    p | omega_matter_now_;
    p | omega_dark_matter_now_;
    p | omega_lambda_now_;
    p | comoving_box_size_;
    p | max_expansion_rate_;
    p | initial_redshift_;
    p | final_redshift_;
  };

  enzo_float hubble_constant_now()   { return hubble_constant_now_; }
  enzo_float omega_matter_now()      { return omega_matter_now_; }
  enzo_float omega_dark_matter_now() { return omega_dark_matter_now_; }
  enzo_float omega_lambda_now()      { return omega_lambda_now_; }
  enzo_float comoving_box_size()     { return comoving_box_size_; }
  enzo_float max_expansion_rate()    { return max_expansion_rate_; }
  enzo_float initial_redshift()      { return initial_redshift_; }
  enzo_float final_redshift()        { return final_redshift_; }
  
  enzo_float initial_time_in_code_units() const
  { return time_from_redshift_ (initial_redshift_); }

  void compute_expansion_factor (enzo_float *a, enzo_float *dadt, enzo_float time) const;
  void compute_expansion_timestep (enzo_float *dt_expansion, enzo_float time) const;

  void get_units
  (enzo_float * density_units,
   enzo_float * length_units,
   enzo_float * temperature_units,
   enzo_float * time_units,
   enzo_float * velocity_units,
   enzo_float time) const;

public: // virtual methods

  virtual std::string type() const { return "cosmology"; }

private: // methods

  enzo_float time_from_redshift_ (enzo_float redshift) const;

private: // attributes
  
  // NOTE: change pup() function whenever attributes change

  enzo_float hubble_constant_now_;
  enzo_float omega_matter_now_;
  enzo_float omega_dark_matter_now_;
  enzo_float omega_lambda_now_;
  enzo_float comoving_box_size_;
  enzo_float max_expansion_rate_;
  enzo_float initial_redshift_;
  enzo_float final_redshift_;

};

#endif /* ENZO_ENZO_PHYSICS_COSMOLOGY_HPP */

