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

  //--------------------------------------------------

public: // virtual functions

  /// Solve the linear system Ax = b
  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw();
  
  /// Return the name of this solver
  virtual std::string name () const
  { return "cg"; }

  //--------------------------------------------------
  
public: // virtual functions

  /// Continuation after global reduction
  template <class T>
  void shift_1(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  template <class T>
  void loop_2a(EnzoBlock * enzo_block) throw();

    /// Continuation after global reduction
  template <class T>
  void loop_2b(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  template <class T>
  void loop_4(EnzoBlock * enzo_block) throw();

  // /// Continuation after global reduction
  // template <class T>
  // void shift_2(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  template <class T>
  void loop_6(EnzoBlock * enzo_block) throw();

  template <class T>
  void end (EnzoBlock * enzo_block, int retval) throw();

  /// Set rz_ by EnzoBlock after reduction
  void set_rz(long double rz) throw()    {  rz_ = rz; }

  /// Set rr_ by EnzoBlock after reduction
  void set_rr(long double rr) throw()    {  rr_ = rr; }

  /// Set rr_new_ by EnzoBlock after reduction
  void set_rz2(long double rz2) throw()  {  rz2_ = rz2; }

  /// Set dy_ by EnzoBlock after reduction
  void set_dy(long double dy) throw()         { dy_ = dy; }

  /// Set bs_ (B sum) by EnzoBlock after reduction
  void set_bs(long double bs) throw()    { bs_ = bs;  }
  /// Set rs_ (R sum) by EnzoBlock after reduction
  void set_rs(long double rs) throw()    { rs_ = rs;  }
  /// Set xs_ (X sum) by EnzoBlock after reduction
  void set_xs(long double xs) throw()    { xs_ = xs;  }

  /// Set bc_ (B count) by EnzoBlock after reduction
  void set_bc(long double bc) throw()    { bc_ = bc;  }

  /// Set iter_ by EnzoBlock after reduction
  void set_iter(int iter) throw()        { iter_ = iter; }

protected: // methods

  template <class T>
  void compute_ (EnzoBlock * enzo_block) throw();

  template <class T>
  void begin_1_() throw();

  void exit_() throw();

  /// Compute local sum of vector elements X_i
  template <class T>
  long double sum_ (const T * X) const throw();

  /// scale the vector by the given scalar Y = a*X
  template <class T>
  void scale_ (T * Y, T a, const T * X) const throw();

  /// return the number of elements of the vector X
  int count_ () const throw();
  
  /// Shift the vector X by a scalar multiple of Y
  /// NOTE includes ghost zones since performed after ghost refresh
  template <class T>
  void shift_ (T * X, const T a, const T * Y) const throw();

protected: // attributes

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

