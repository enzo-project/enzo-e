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

  EnzoSolverCg (std::string name,
		std::string field_x,
		std::string field_b,
		int monitor_iter,
		int restart_cycle,
		int solve_type,
		int min_level,
		int max_level,
		int iter_max, 
		double res_tol,
		int index_precon);

  /// Constructor
  EnzoSolverCg() throw()
  : Solver(), 
    A_(NULL),
    index_precon_(-1),
    iter_max_(0), 
    ir_(-1), id_(-1), iy_(-1), iz_(-1),
    nx_(0),ny_(0),nz_(0),
    mx_(0),my_(0),mz_(0),
    gx_(0),gy_(0),gz_(0),
    iter_(0),
    res_tol_(0.0),
    rr0_(0),
    rr_min_(0),rr_max_(0),
    rr_(0.0), rz_(0.0), rz2_(0.0), dy_(0.0), bs_(0.0), rs_(0.0), xs_(0.0),
    bc_(0.0),
    local_(false)
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverCg);

  /// Charm++ PUP::able migration constructor
  EnzoSolverCg (CkMigrateMessage *m)
    : Solver(m), 
      A_(NULL),
      index_precon_(-1),
      iter_max_(0), 
      ir_(-1), id_(-1), iy_(-1), iz_(-1),
      nx_(0),ny_(0),nz_(0),
      mx_(0),my_(0),mz_(0),
      gx_(0),gy_(0),gz_(0),
      iter_(0),
      res_tol_(0.0),
      rr0_(0),
      rr_min_(0),rr_max_(0),
      rr_(0.0), rz_(0.0), rz2_(0.0), dy_(0.0), bs_(0.0), rs_(0.0), xs_(0.0),
      bc_(0.0),
      local_(false)
  {}

  /// Assignment operator
  EnzoSolverCg & operator= (const EnzoSolverCg & EnzoSolverCg) throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  //--------------------------------------------------

public: // virtual functions

  /// Solve the linear system Ax = b
  virtual void apply ( std::shared_ptr<Matrix> A, Block * block) throw();
  
  /// Type of this solver
  virtual std::string type() const { return "cg"; }

  //--------------------------------------------------
  
public: // virtual functions

  /// Continuation after global reduction
  void shift_1(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  void loop_0a(EnzoBlock * enzo_block, CkReductionMsg *) throw();

    /// Continuation after global reduction
  void loop_0b(EnzoBlock * enzo_block, CkReductionMsg *) throw();

  /// Continuation after global reduction
  void loop_2a(EnzoBlock * enzo_block) throw();

    /// Continuation after global reduction
  void loop_2b(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  void loop_4(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  void loop_6(EnzoBlock * enzo_block) throw();

  void end (EnzoBlock * enzo_block, int retval) throw();

  /// Set rz_ by EnzoBlock after reduction
  void set_rz(double rz) throw()    {  rz_ = rz; }

  /// Set rr_ by EnzoBlock after reduction
  void set_rr(double rr) throw()    {  rr_ = rr; }

  /// Set rr_new_ by EnzoBlock after reduction
  void set_rz2(double rz2) throw()  {  rz2_ = rz2; }

  /// Set dy_ by EnzoBlock after reduction
  void set_dy(double dy) throw()         { dy_ = dy; }

  /// Set bs_ (B sum) by EnzoBlock after reduction
  void set_bs(double bs) throw()    { bs_ = bs;  }
  /// Set rs_ (R sum) by EnzoBlock after reduction
  void set_rs(double rs) throw()    { rs_ = rs;  }
  /// Set xs_ (X sum) by EnzoBlock after reduction
  void set_xs(double xs) throw()    { xs_ = xs;  }

  /// Set bc_ (B count) by EnzoBlock after reduction
  void set_bc(double bc) throw()    { bc_ = bc;  }

  /// Set iter_ by EnzoBlock after reduction
  void set_iter(int iter) throw()        { iter_ = iter; }

protected: // methods

  void compute_ (EnzoBlock * enzo_block) throw();

  void begin_1_() throw();

  /// Allocate temporary Fields
  void allocate_temporary_(Field field, Block * block = NULL)
  {
    field.allocate_temporary(id_);
    field.allocate_temporary(ir_);
    field.allocate_temporary(iy_);
    field.allocate_temporary(iz_);
  }

  /// Dellocate temporary Fields
  void deallocate_temporary_(Field field, Block * block = NULL)
  {
    field.deallocate_temporary(id_);
    field.deallocate_temporary(ir_);
    field.deallocate_temporary(iy_);
    field.deallocate_temporary(iz_);
  }

  /// Serial CG solver if local_ == true
  void local_cg_ (EnzoBlock * enzo_block);

  /// Apply boundary conditions for the Field on the local block
  void refresh_local_(int ix, EnzoBlock * enzo_block);

  /// Shift Field so that sum(x) == 0
  void shift_local_(int ix, EnzoBlock * enzo_block);
  
  void monitor_output_(EnzoBlock *);
  
protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Matrix
  std::shared_ptr<Matrix> A_;

  /// Preconditioner (-1 if none)
  int index_precon_;

  /// Maximum number of Cg iterations
  int iter_max_;

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

  /// Convergence tolerance on the residual reduction rz_ / rz0_
  double res_tol_;

  /// Initial residual
  double rr0_;

  /// Minimum residual
  double rr_min_;

  /// Maximum residual
  double rr_max_;

  /// dot (R_i,R_i)
  double rr_;

  /// dot (R_i,Z_i)
  double rz_;

  /// dot (R_i+1,Z_i+1)
  double rz2_;

  /// dot (D,Y)
  double dy_;

  /// sum of elements B(i) for singular systems
  double bs_;
  /// sum of elements R(i) for singular systems
  double rs_;
  /// sum of elements X(i) for singular systems
  double xs_;

  /// count of elements B(i) for singular systems
  double bc_;

  /// Whether to solve on a standalone Block, e.g. for MG coarse solver
  bool local_;
};

#endif /* ENZO_ENZO_SOLVER_CG_HPP */

