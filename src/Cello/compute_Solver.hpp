// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Solver.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Compute] Declaration for the Solver class

#ifndef COMPUTE_SOLVER_HPP
#define COMPUTE_SOLVER_HPP

class Solver : public PUP::able 
{
  /// @class    Solver
  /// @ingroup  Compute
  /// @brief    [\ref Solver] Interface to a linear solver

public: // interface

  /// Create a new Solver
  Solver () throw()
  {}

  /// Destructor
  virtual ~Solver() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_abstract(Solver);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);
  }

public: // virtual functions

  /// Solve the linear system Ax = b

  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw() = 0; 

protected: // functions

};

#endif /* COMPUTE_SOLVER_HPP */
