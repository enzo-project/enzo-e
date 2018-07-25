// See LICENSE_CELLO file for license and copyright information

/// @file     compute_SolverNull.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2018-02-26
/// @brief    [\ref Compute] Declaration for the SolverNull class

#ifndef COMPUTE_SOLVER_NULL_HPP
#define COMPUTE_SOLVER_NULL_HPP

#include <cstring>

class SolverNull : public Solver
{
  /// @class    SolverNull
  /// @ingroup  Compute
  /// @brief    [\ref SolverNull] Interface to a linear SolverNull

public: // interface

  /// Create a new SolverNull
  SolverNull (std::string name,
	      int monitor_iter,
	      int restart_cycle,
	      int min_level = 0,
	      int max_level = std::numeric_limits<int>::max()) throw()
    : Solver(name,monitor_iter,restart_cycle,min_level,max_level)
  {}

  /// Create an uninitialized SolverNull
  SolverNull () throw()
  : Solver()
  {}

  /// Destructor
  virtual ~SolverNull() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_decl(SolverNull);

  SolverNull (CkMigrateMessage *m)
    : Solver(m)
  { }
  
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual functions

  /// Solve the linear system Ax = b
  virtual void apply ( std::shared_ptr<Matrix> A,
		       int ix, int ib, Block * block) throw();

  /// Return the name of this SolverNull
  virtual std::string name () const
  { return "null"; }

protected: // functions
protected: // attributes

};

#endif /* COMPUTE_SOLVER_NULL_HPP */
