// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverCg.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-11-08
/// @brief    [\ref Enzo] Declaration of the EnzoSolverCg class

#ifndef ENZO_ENZO_SOLVER_CG_HPP
#define ENZO_ENZO_SOLVER_CG_HPP

class EnzoSolverCg : public Solver {

  /// @class    EnzoSolverCg
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  EnzoSolverCg (const FieldDescr * field_descr,
		Matrix * A,
		int rank,
		int iter_max, 
		double res_tol,
		int monitor_iter,
		bool is_singular,
		bool diag_precon);

  /// Constructor
  EnzoSolverCg() throw()
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverCg);

  /// Charm++ PUP::able migration constructor
  EnzoSolverCg (CkMigrateMessage *m)
  {}

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

  /// Matrix
  Matrix * A_;

  /// Solution and right-hand-side fields
  int ix_;
  int ib_;

  /// Preconditioner
  Matrix * M_;

  /// Whether you need to subtract of the nullspace of A from b, e.g. fully
  /// periodic or Neumann problems
  bool is_singular_;

  /// Dimensionality of the problem
  int rank_;

  /// Maximum number of Cg iterations
  int iter_max_;

  /// Convergence tolerance on the residual reduction rz_ / rz0_
  double res_tol_;

  /// How often to display progress
  int monitor_iter_;

  /// Initial residual
  long double rr0_;

  /// Minimum residual
  long double rr_min_;

  /// Maximum residual
  long double rr_max_;

  /// Density and potential field id's

  int idensity_;
  int ipotential_;

  /// CG vector id's
  int ir_;
  int id_;
  int iy_;
  int iz_;

  /// Block field attributes
  int nx_,ny_,nz_;
  int mx_,my_,mz_;
  int gx_,gy_,gz_;

  /// Current CG iteration
  int iter_;

  /// dot (R_i,R_i)
  long double rr_;

  /// dot (R_i,Z_i)
  long double rz_;

  /// dot (R_i+1,Z_i+1)
  long double rz2_;

  /// dot (D,Y)
  long double dy_;

  /// sum of elements B(i) for singular systems
  long double bs_;
  /// sum of elements R(i) for singular systems
  long double rs_;
  /// sum of elements X(i) for singular systems
  long double xs_;

  /// count of elements B(i) for singular systems
  long double bc_;

  /// matvec refresh index
  int id_refresh_matvec_;

};

#endif /* ENZO_ENZO_SOLVER_CG_HPP */

