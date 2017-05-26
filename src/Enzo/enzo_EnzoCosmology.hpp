// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCosmology.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-05-23
/// @brief    [\ref Enzo] Declaration of the EnzoCosmology class

#ifndef ENZO_ENZO_COSMOLOGY_HPP
#define ENZO_ENZO_COSMOLOGY_HPP

class EnzoCosmology {

  /// @class    EnzoCosmology
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoCosmology() throw()
  : hubble_constant_now_(0.701),
    omega_matter_now_(0.279),
    omega_dark_matter_now_(-1.0),
    omega_lambda_now_(0.721),
    comoving_box_size_(64),
    max_expansion_rate_(0.01),
    initial_redshift_(20.0),
    final_redshift_(0.0)
  {  }

  /// Constructor
  EnzoCosmology
  (
   double hubble_constant_now,
   double omega_matter_now,
   double omega_dark_matter_now,
   double omega_lambda_now,
   double comoving_box_size,
   double max_expansion_rate,
   double initial_redshift,
   double final_redshift
   )
    : hubble_constant_now_(hubble_constant_now),
      omega_matter_now_(omega_matter_now),
      omega_dark_matter_now_(omega_dark_matter_now),
      omega_lambda_now_(omega_lambda_now),
      comoving_box_size_(comoving_box_size),
      max_expansion_rate_(max_expansion_rate),
      initial_redshift_(initial_redshift),
      final_redshift_(final_redshift)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    
    p | hubble_constant_now_;
    p | omega_matter_now_;
    p | omega_dark_matter_now_;
    p | omega_lambda_now_;
    p | comoving_box_size_;
    p | max_expansion_rate_;
    p | initial_redshift_;
    p | final_redshift_;
  };

  double hubble_constant_now()
  { return hubble_constant_now_; }
  double omega_matter_now()
  { return omega_matter_now_; }
  double omega_dark_matter_now()
  { return omega_dark_matter_now_; }
  double omega_lambda_now()
  { return omega_lambda_now_; }
  double comoving_box_size()
  { return comoving_box_size_; }
  double max_expansion_rate()
  { return max_expansion_rate_; }
  double initial_redshift()
  { return initial_redshift_; }
  double final_redshift()
  { return final_redshift_; }
  
  double initial_time_in_code_units()
  { return time_from_redshift_ (initial_redshift_); }

private: // methods

  double time_from_redshift_ (double redshift) const;

private: // attributes
  
  // NOTE: change pup() function whenever attributes change

  double hubble_constant_now_;
  double omega_matter_now_;
  double omega_dark_matter_now_;
  double omega_lambda_now_;
  double comoving_box_size_;
  double max_expansion_rate_;
  double initial_redshift_;
  double final_redshift_;

};

#endif /* ENZO_ENZO_COSMOLOGY_HPP */

