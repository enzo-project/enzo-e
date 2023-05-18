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
    omega_baryon_now_(0.279),
    omega_cdm_now_(0.0),
    omega_lambda_now_(0.721),
    comoving_box_size_(64),
    max_expansion_rate_(0.01),
    initial_redshift_(20.0),
    final_redshift_(0.0),
    cosmo_a_(0.0),
    cosmo_dadt_(0.0),
    current_redshift_(-1.0)
  {
  }

  /// Constructor
  EnzoPhysicsCosmology
  (
   enzo_float hubble_constant_now,
   enzo_float omega_matter_now,
   enzo_float omega_baryon_now,
   enzo_float omega_cdm_now,
   enzo_float omega_lambda_now,
   enzo_float comoving_box_size,
   enzo_float max_expansion_rate,
   enzo_float initial_redshift,
   enzo_float final_redshift
   )
    : Physics(),
      hubble_constant_now_(hubble_constant_now),
      omega_matter_now_(omega_matter_now),
      omega_baryon_now_(omega_baryon_now),
      omega_cdm_now_(omega_cdm_now),
      omega_lambda_now_(omega_lambda_now),
      comoving_box_size_(comoving_box_size),
      max_expansion_rate_(max_expansion_rate),
      initial_redshift_(initial_redshift),
      final_redshift_(final_redshift),
      cosmo_a_(0.0),
      cosmo_dadt_(0.0),
      current_redshift_(-1.0)
  {
    ASSERT3 ("EnzoPhysicsCosmology::EnzoPhysicsCosmology()",
	     "omega_matter_now (%g) must equal "
	     "omega_cdm_now (%g) + omega_baryon_now (%g)",
	     omega_matter_now_,omega_cdm_now_,omega_baryon_now_,
	     std::abs
	     (omega_matter_now_-(omega_cdm_now_+omega_baryon_now_)) < 1e-7);
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoPhysicsCosmology);

  /// CHARM++ migration constructor
  EnzoPhysicsCosmology(CkMigrateMessage *m)
    : Physics (m),
      hubble_constant_now_(0.701),
      omega_matter_now_(0.279),
      omega_baryon_now_(1.0),
      omega_cdm_now_(0.0),
      omega_lambda_now_(0.721),
      comoving_box_size_(64),
      max_expansion_rate_(0.01),
      initial_redshift_(20.0),
      final_redshift_(0.0),
      cosmo_a_(0.0),
      cosmo_dadt_(0.0),
      current_redshift_(-1.0)
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
    p | omega_baryon_now_;
    p | omega_cdm_now_;
    p | omega_lambda_now_;
    p | comoving_box_size_;
    p | max_expansion_rate_;
    p | initial_redshift_;
    p | final_redshift_;

    p | cosmo_a_;
    p | cosmo_dadt_;
    p | current_redshift_;

  };

  enzo_float hubble_constant_now()   { return hubble_constant_now_; }
  enzo_float omega_matter_now()      { return omega_matter_now_; }
  enzo_float omega_baryon_now()      { return omega_baryon_now_; }
  enzo_float omega_cdm_now()         { return omega_cdm_now_; }
  enzo_float omega_lambda_now()      { return omega_lambda_now_; }
  enzo_float comoving_box_size()     { return comoving_box_size_; }
  enzo_float max_expansion_rate()    { return max_expansion_rate_; }
  enzo_float current_redshift()      { return current_redshift_; }
  enzo_float initial_redshift()      { return initial_redshift_; }
  enzo_float final_redshift()        { return final_redshift_; }

  /// Set Cosmology parameters (used by testing only)
  void set_hubble_constant_now(enzo_float value)
  { hubble_constant_now_=value; }
  void set_omega_matter_now(enzo_float value)
  { omega_matter_now_=value; }
  void set_omega_baryon_now(enzo_float value)
  { omega_baryon_now_=value; }
  void set_omega_cdm_now(enzo_float value)
  { omega_cdm_now_=value; }
  void set_omega_lambda_now(enzo_float value)
  { omega_lambda_now_=value; }
  void set_comoving_box_size(enzo_float value)
  { comoving_box_size_=value; }
  void set_max_expansion_rate(enzo_float value)
  { max_expansion_rate_=value; }
  void set_initial_redshift(enzo_float value)
  { initial_redshift_=value; }
  void set_final_redshift(enzo_float value)
  { final_redshift_=value; }

  enzo_float initial_time_in_code_units() const
  { return time_from_redshift (initial_redshift_); }
  enzo_float time_from_redshift (enzo_float redshift) const;
  enzo_float redshift_from_time(enzo_float time) const
  {
    enzo_float cosmo_a,cosmo_dadt;
    compute_expansion_factor (&cosmo_a,&cosmo_dadt,time);
    return (1.0 + initial_redshift_) / cosmo_a - 1.0;
  }

  void set_current_time (enzo_float time)
  {
    update_expansion_factor (time);
    current_redshift_ = (1 + initial_redshift_)/cosmo_a_ - 1;

  }

  void set_current_redshift (enzo_float redshift)
  {
    set_current_time(time_from_redshift(redshift));
  }

  void update_expansion_factor(enzo_float time)
  {
    compute_expansion_factor(&cosmo_a_,&cosmo_dadt_,time);
  }

  void compute_expansion_timestep
  (enzo_float *dt_expansion, enzo_float time) const;

  void compute_expansion_factor
  (enzo_float *cosmo_a, enzo_float *cosmo_dadt, enzo_float time) const;

  /// Return density unit at current time / redshift in terms of `g / cm^3`
  /// (requires set_current_time()).
  double density_units() const
  {

    /* Density unit is defined so that the mean physical matter density of the universe is 1.
     This can be written as:

     `rho_bar_m_0 * (1 + z)^3 / rho_unit = 1`

     Where `rho_bar_m_0` is the mean physical matter density at redshift 0, `z` is redshift,
     and `rho_unit` is the density unit in `g / cm^3`.

     We can then write this as:

     `rho_unit = `Omega_m_0 * rho_crit_0 * (1 + z)^3`,

     where `Omega_m_0` is a cosmological parameter which can be set in the input parameter file,
     and `rho_crit_0` is defined as `3 * H_0^2 / (8 * pi * G)`, where `H_0` is the expansion rate
     of the universe at redshift 0 (with dimensions of inverse time), and `G` is the
     gravitational constant. `H_0` can be written as `H_0 / h * h`, where `h` is given the name
     `hubble_constant_now_` in the code below. Putting this all together gives the following
     expression:

     `rho_unit = (3 / 8 * pi) * (H_0 / h)^2 / G * Omega_m_0 * h^2 * (1+z)^3`.

     Note: the previous version of this function was equivalent to:

     `rho_unit = 1.8788e-29 * Omega_m_0 * h^2 * (1+z)^3`.

     but (after plugging in values), this version is equivalent to:

     `rho_unit = 1.8784710838431666e-29 * Omega_m_0 * h^2 * (1+z)^3`.

     The reason for this inconsistency is unclear. It is possibly due to inconsistent values
     for `G`.

   */

    return (3.0 / (8.0 * cello::pi)) * enzo_constants::H0_over_h * enzo_constants::H0_over_h /
      enzo_constants::grav_constant * omega_matter_now_ *
      hubble_constant_now_ * hubble_constant_now_ *
      (1 + current_redshift_) * (1 + current_redshift_) * (1 + current_redshift_);
  }

  /// Return current mass unit in terms of grams (requires set_current_time())
  double mass_units() const
  {
    return density_units() * length_units() * length_units() * length_units();
  }

  /// Return current length unit in terms of cm (requires set_current_time())
  double length_units() const
  {
    return enzo_constants::Mpc_cm * comoving_box_size_ / hubble_constant_now_ /
      (1.0 + current_redshift_);
  }

  /// Return time unit in terms of seconds (requires set_current_time())
  double time_units() const
  {

    /* Time unit is defined so that the following is true (which has the effect of simplifying
      Poisson's equation):

      `4 * pi * G * rho_bar_m_com * time_unit^2 / a_unit^3 = 1`,

      Where `G` is the gravitational constant.`rho_bar_m_com` is the mean comoving matter density
      of the universe, `time_unit` is the time unit in seconds, and `a_unit` is the "unit 
      cosmological scale factor". To understand where the factor of `a_unit^3` come from, see
      Equations 8, 17, 18 of Greg L. Bryan et al 2014 ApJS 211 19, and note that Laplacian(`phi`) 
      has dimensions of `time^(-2) * a^(2).

      In Enzo-E, `a_unit` is defined so that `a` is 1 at the initial redshift (`z_i`), so that:
      `(1+z_i)^(-1) / a_unit = 1`, which means that `a_unit = 1 / (1 + z_i)`.

      `rho_bar_m_0` can be written as `Omega_m_0 * rho_crit_0`, where `Omega_m_0` is a 
      cosmological parameter which can be set in the input parameter file,
      and `rho_crit_0` is defined as `3 * H_0^2 / (8 * pi * G)`, where `H_0` is the expansion rate
      of the universe at redshift 0 (with dimensions of inverse time). Plugging this all in gives:

      `3 / 2 * Omega_m_0 * H_0^2 * time_unit^2 * (1 + z_i)^3 = 1`.

      Equivalently, this is the "free-fall time" at `z = z_i`.

      After some rearrangement, we get the following expression:

      `time_unit = sqrt(2/3) * h / H0 / (h * sqrt(Omega_m_0 * (1 + z)^3))`/

   */

    return sqrt(2.0/3.0) / (enzo_constants::H0_over_h * hubble_constant_now_ *
			    sqrt(omega_matter_now_ *
				 (1 + initial_redshift_) *
				 (1 + initial_redshift_) *
				 (1 + initial_redshift_)));
  }

  /// Return velocity unit in cm/s
  double velocity_units() const
  {
    /* Equation 17 of Greg L. Bryan et al 2014 ApJS 211 19, shows that the dimensions of velocity
       are `a` times `comoving length` over `time`.
       As such, this expression is equivalent to:

       `velocity_unit` = `a_unit` * `length_unit / ((1 + z) * `time_unit`)

       where `z` is the current redshift, and `a_unit` is (1 + `z_i`), with `z_i` being the 
       initial redshift.

       The `1.0e7` term is the result of multiplying `Mpc_cm` by `H0_over_h`.
    */
    return 1.0e7 *  sqrt(3.0/2.0) * comoving_box_size_ * sqrt(omega_matter_now_) *
      sqrt(1 + initial_redshift_);
  }

  void print () const
  {
    CkPrintf ("EnzoPhysicsCosmology hubble_constant_now = %g\n",hubble_constant_now_);
    CkPrintf ("EnzoPhysicsCosmology omega_matter_now    = %g\n",omega_matter_now_);
    CkPrintf ("EnzoPhysicsCosmology omega_baryon_now    = %g\n",omega_baryon_now_);
    CkPrintf ("EnzoPhysicsCosmology omega_cdm_now       = %g\n",omega_cdm_now_);
    CkPrintf ("EnzoPhysicsCosmology omega_lambda_now    = %g\n",omega_lambda_now_);
    CkPrintf ("EnzoPhysicsCosmology comoving_box_size   = %g\n",comoving_box_size_);
    CkPrintf ("EnzoPhysicsCosmology max_expansion_rate  = %g\n",max_expansion_rate_);
    CkPrintf ("EnzoPhysicsCosmology initial_redshift    = %g\n",initial_redshift_);
    CkPrintf ("EnzoPhysicsCosmology final_redshift      = %g\n",final_redshift_);
    CkPrintf ("EnzoPhysicsCosmology cosmo_a             = %g\n",cosmo_a_);
    CkPrintf ("EnzoPhysicsCosmology cosmo_dadt          = %g\n",cosmo_dadt_);
    CkPrintf ("EnzoPhysicsCosmology current_redshift    = %g\n",current_redshift_);
    fflush(stdout);
  }

public: // virtual methods

  virtual std::string type() const { return "cosmology"; }

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  // Constant parameters
  enzo_float hubble_constant_now_;
  enzo_float omega_matter_now_;
  enzo_float omega_baryon_now_;
  enzo_float omega_cdm_now_;
  enzo_float omega_lambda_now_;
  enzo_float comoving_box_size_;
  enzo_float max_expansion_rate_;
  enzo_float initial_redshift_;
  enzo_float final_redshift_;

  // Time-dependent parameters
  enzo_float cosmo_a_;
  enzo_float cosmo_dadt_;
  enzo_float current_redshift_;

};

#endif /* ENZO_ENZO_PHYSICS_COSMOLOGY_HPP */
