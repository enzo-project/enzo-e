// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodNull.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Problem] Declaration of MethodNull "null" solver

#ifndef PROBLEM_METHOD_NULL_HPP
#define PROBLEM_METHOD_NULL_HPP

class MethodNull : public Method {

  /// @class    MethodNull
  /// @ingroup  Problem
  ///
  /// @brief [\ref Problem] Null method, used for forcing refresh
  /// before other methods at the start of each cycle

public: // interface

  /// Create a new MethodNull object
  MethodNull ( double dt )
    : Method(), dt_(dt)
  {
    init_refresh_();
  }

  MethodNull()
    : Method(), dt_ (std::numeric_limits<double>::max())
  {
    init_refresh_();
  }

  /// Initialize refresh
  void init_refresh_();
  
  /// Charm++ PUP::able declarations
  PUPable_decl(MethodNull);
  
  /// Charm++ PUP::able migration constructor
  MethodNull (CkMigrateMessage *m)
    : Method (m), dt_(0.0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; Method::pup(p); p | dt_; }
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "null"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw()
  { return dt_; }

protected: // attributes

  /// Time step
  double dt_;
};

#endif /* PROBLEM_METHOD_NULL_HPP */
