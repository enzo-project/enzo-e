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
  EnzoSolverDiagonal(std::string name,
		     std::string field_x,
		     std::string field_b,
		     int monitor_iter,
		     int restart_cycle,
		     int solve_type) throw()
    : Solver(name,field_x,field_b,monitor_iter,restart_cycle,solve_type)
  { }

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
  virtual void apply ( std::shared_ptr<Matrix> A, Block * block) throw();
  
  /// Type of this solver
  virtual std::string type() const { return "diagonal"; }

  //--------------------------------------------------
  
public: // virtual functions

protected: // methods

  void compute_ (std::shared_ptr<Matrix> A, Block * block) throw();

protected: // attributes

};

#endif /* ENZO_ENZO_SOLVER_DIAGONAL_HPP */

