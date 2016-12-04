// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverCg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-08
/// @brief    [\ref Enzo] Declaration of the EnzoSolverCg class

#ifndef ENZO_ENZO_SOLVER_CG_HPP
#define ENZO_ENZO_SOLVER_CG_HPP

class EnzoSolverCg : Solver {

  /// @class    EnzoSolverCg
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoSolverCg() throw();

  /// Copy constructor
  EnzoSolverCg(const EnzoSolverCg & EnzoSolverCg) throw();

  /// Assignment operator
  EnzoSolverCg & operator= (const EnzoSolverCg & EnzoSolverCg) throw();

  /// Destructor
  ~EnzoSolverCg() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual functions

  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw();
  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* ENZO_ENZO_SOLVER_CG_HPP */

