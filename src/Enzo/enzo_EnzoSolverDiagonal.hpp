// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverDiagonal.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-03-01
/// @brief    [\ref Enzo] Declaration of the EnzoSolverDiagonal class

#ifndef ENZO_ENZO_SOLVER_DIAGONAL_HPP
#define ENZO_ENZO_SOLVER_DIAGONAL_HPP

class EnzoSolverDiagonal : public Solver {

  /// @class    EnzoSolverDiagonal
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoSolverDiagonal() throw()
  : Solver()
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverDiagonal);

  /// Charm++ PUP::able migration constructor
  EnzoSolverDiagonal (CkMigrateMessage *m)
    : Solver(m)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    Solver::pup(p);
  };

  //--------------------------------------------------

public: // virtual functions

  /// Solve the linear system Ax = b
  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw();
  
  /// Return the name of this solver
  virtual std::string name () const
  { return "diagonal"; }

  //--------------------------------------------------
  
public: // virtual functions

protected: // methods

  template <class T>
  void compute_ (Matrix * A, int ix, int ib, Block * block) throw();

protected: // attributes

};

#endif /* ENZO_ENZO_SOLVER_DIAGONAL_HPP */

