// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverJacobi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-04-30 18:45:58
/// @brief    [\ref Enzo] Declaration of the EnzoSolverJacobi class

#ifndef ENZO_ENZO_SOLVER_JACOBI_HPP
#define ENZO_ENZO_SOLVER_JACOBI_HPP

class EnzoSolverJacobi : public Solver {

  /// @class    EnzoSolverJacobi
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoSolverJacobi(const FieldDescr * field_descr,
		   double weight=1.0,
		   int iter_max = 1) throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverJacobi);

  /// Charm++ PUP::able migration constructor
  EnzoSolverJacobi (CkMigrateMessage *m)
    : A_(NULL),
      ix_(0),
      ib_(0),
      ir_(0),
      id_(0),
      w_(0),
      n_(0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { 
    TRACEPUP;
    Solver::pup(p);

    p | A_;
    p | ix_;
    p | ib_;
    p | ir_;
    p | id_;
    p | w_;
    p | n_;
  }

public: // virtual functions

  /// Solve the linear system Ax = b
  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw();

  /// Return the name of this solver
  virtual std::string name () const
  { return "jacobi"; }

public: // functions

  /// Continue after refresh to perform Jacobi update
  void compute (Block * block);
  
protected: // functions

  /// Implementation of solver() for given precision 
  template <typename T>
  void apply_(Block * block);

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Matrix A for smoothing A*X = B
  Matrix * A_;
  
  /// Field index for vector X
  int ix_;

  /// Field index for rhs B
  int ib_;

  /// Field index for residual R
  int ir_;

  /// Field index for matrix diagonal D
  int id_;

  /// Weighting
  double w_;

  /// Number of iterations
  int n_;

};

#endif /* ENZO_ENZO_SOLVER_JACOBI_HPP */

