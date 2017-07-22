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
    omega_lambda_now_(0.721),
    comoving_box_size_(64),
    max_expansion_rate_(0.01),
    initial_redshift_(20.0),
    final_redshift_(0.0),
    a_(0.0),
    dadt_(0.0),
    current_redshift_(-1.0)
  {  }

  /// Constructor
  EnzoPhysicsCosmology
  (
   enzo_float hubble_constant_now,
   enzo_float omega_matter_now,
   enzo_float omega_lambda_now,
   enzo_float comoving_box_size,
   enzo_float max_expansion_rate,
   enzo_float initial_redshift,
   enzo_float final_redshift
   )
    : Physics(),
      hubble_constant_now_(hubble_constant_now),
      omega_matter_now_(omega_matter_now),
      omega_lambda_now_(omega_lambda_now),
      comoving_box_size_(comoving_box_size),
      max_expansion_rate_(max_expansion_rate),
      initial_redshift_(initial_redshift),
      final_redshift_(final_redshift),
      a_(0.0),
      dadt_(0.0),
      current_redshift_(-1.0)
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoPhysicsCosmology);

  /// CHARM++ migration constructor
  EnzoPhysicsCosmology(CkMigrateMessage *m)
    : Physics (m)
  {}

  /// Virtual destructor
  virtual ~EnzoPhysicsCosmology()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    Physics::pup(p);
    
    p | hubble_constant_now_;
    p | omega_matter_now_;
    p | omega_lambda_now_;
    p | comoving_box_size_;
    p | max_expansion_rate_;
    p | initial_redshift_;
    p | final_redshift_;

    p | a_;
    p | dadt_;
    p | current_redshift_;
  };

  enzo_float hubble_constant_now()   { return hubble_constant_now_; }
  enzo_float omega_matter_now()      { return omega_matter_now_; }
  enzo_float omega_lambda_now()      { return omega_lambda_now_; }
  enzo_float comoving_box_size()     { return comoving_box_size_; }
  enzo_float max_expansion_rate()    { return max_expansion_rate_; }
  enzo_float initial_redshift()      { return initial_redshift_; }
  enzo_float final_redshift()        { return final_redshift_; }

  enzo_float set_hubble_constant_now(enzo_float value)
  { hubble_constant_now_=value; }
  enzo_float set_omega_matter_now(enzo_float value)
  { omega_matter_now_=value; }
  enzo_float set_omega_lambda_now(enzo_float value)
  { omega_lambda_now_=value; }
  enzo_float set_comoving_box_size(enzo_float value)
  { comoving_box_size_=value; }
  enzo_float set_max_expansion_rate(enzo_float value)
  { max_expansion_rate_=value; }
  enzo_float set_initial_redshift(enzo_float value)
  { initial_redshift_=value; }
  enzo_float set_final_redshift(enzo_float value)
  { final_redshift_=value; }
  
  enzo_float initial_time_in_code_units() const
  { return time_from_redshift (initial_redshift_); }
  enzo_float time_from_redshift (enzo_float redshift) const;
  enzo_float redshift_from_time(enzo_float time) const
  {
    enzo_float a,dadt;
    compute_expansion_factor (&a,&dadt,time);
    return (1.0 + initial_redshift_) / a - 1.0;
  }

  void set_current_time (enzo_float time)
  {
    update_expansion_factor (time);
    current_redshift_ = (1 + initial_redshift_)/a_ - 1;
    
  }
  void update_expansion_factor(enzo_float time)
  {
    compute_expansion_factor(&a_,&dadt_,time);
  }
      
  void compute_expansion_timestep
  (enzo_float *dt_expansion, enzo_float time) const;

  void compute_expansion_factor
  (enzo_float *a, enzo_float *dadt, enzo_float time) const;
  
  /// Return current mass units scaling (requires set_current_time())
  double mass_units() const
  {
    double density = 1.88e-29*omega_matter_now_*
      pow(hubble_constant_now_,2)*
      pow(1 + current_redshift_,3);
    double length = length_units();
    return density * length * length * length;
  }

  /// Return current length units scaling (requires set_current_time())
  double length_units() const
  {
    return 3.086e24*comoving_box_size_/hubble_constant_now_/
      (1.0 + current_redshift_);
  }

  /// Return current time units (requires set_current_time())
  double time_units() const
  {
    return 2.519445e17/sqrt(omega_matter_now_)/hubble_constant_now_/
      pow(1.0 + initial_redshift_,1.5);
  }

  void print () const
  {
    CkPrintf ("hubble_constant_now_   = %g\n",hubble_constant_now_);
    CkPrintf ("omega_matter_now_      = %g\n",omega_matter_now_);
    CkPrintf ("omega_lambda_now_      = %g\n",omega_lambda_now_);
    CkPrintf ("comoving_box_size_     = %g\n",comoving_box_size_);
    CkPrintf ("max_expansion_rate_    = %g\n",max_expansion_rate_);
    CkPrintf ("initial_redshift_      = %g\n",initial_redshift_);
    CkPrintf ("final_redshift_        = %g\n",final_redshift_);
    CkPrintf ("a_                     = %g\n",a_);
    CkPrintf ("dadt_                  = %g\n",dadt_);
    CkPrintf ("current_redshift_      = %g\n",current_redshift_);
  }

public: // virtual methods

  virtual std::string type() const { return "cosmology"; }

protected: // attributes
  
  // NOTE: change pup() function whenever attributes change

  // Constant parameters
  enzo_float hubble_constant_now_;
  enzo_float omega_matter_now_;
  enzo_float omega_lambda_now_;
  enzo_float comoving_box_size_;
  enzo_float max_expansion_rate_;
  enzo_float initial_redshift_;
  enzo_float final_redshift_;

  // Time-dependent parameters
  enzo_float a_;
  enzo_float dadt_;
  enzo_float current_redshift_;

};

#endif /* ENZO_ENZO_PHYSICS_COSMOLOGY_HPP */

