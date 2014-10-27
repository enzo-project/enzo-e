// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Compute.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Compute] Declaration for the Compute class

#ifndef COMPUTE_COMPUTE_HPP
#define COMPUTE_COMPUTE_HPP

class Compute : public PUP::able 
{
  /// @class    Compute
  /// @ingroup  Compute
  /// @brief    [\ref Compute] Interface to an application compute / analysis / visualization function.

public: // interface

  /// Create a new Compute
  Compute () throw()
  {}

  /// Destructor
  virtual ~Compute() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_abstract(Compute);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);
  }

public: // virtual functions

  /// Perform the computation on the CommBlock

  virtual void compute ( CommBlock * comm_block) throw() = 0; 

protected: // functions

};

#endif /* COMPUTE_COMPUTE_HPP */
