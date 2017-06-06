// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Physics.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Problem] Declaration for the Physics component

#ifndef PROBLEM_PHYSICS_HPP
#define PROBLEM_PHYSICS_HPP

class Hierarchy;

class Physics : public PUP::able 
{

  /// @class    Physics
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Encapsulate physics, e.g. cosmological parameters, etc.

public: // interface

  /// empty constructor for charm++ pup()
  Physics() throw()
  {}

  /// CHARM++ PUP::able declaration
  PUPable_decl(Physics);

  /// CHARM++ migration constructor for PUP::able
  Physics (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
  }

  /// Name of this physics type
  virtual std::string type () const
  { return "undefined"; }
};

#endif /* PROBLEM_PHYSICS_HPP */
