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
  EnzoSolverJacobi(std::string name,
		   std::string field_x,
		   std::string field_b,
		   int monitor_iter,
		   int restart_cycle,
		   int solve_type,
		   double weight=1.0,
		   int iter_max = 1) throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverJacobi);

  /// Charm++ PUP::able migration constructor
  EnzoSolverJacobi (CkMigrateMessage *m)
    : Solver(m),
      A_(NULL),
      ir_(-1),
      id_(-1),
      w_(0),
      i_iter_(-1),
      n_(0)
#ifdef NEW_REFRESH
  , ir_smooth_(-1)
#endif    
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { 
    TRACEPUP;
    Solver::pup(p);

    //    p | A_;
    p | ir_;
    p | id_;
    p | w_;
    p | i_iter_;
    p | n_;
#ifdef NEW_REFRESH
    p | ir_smooth_;
#endif    
  }

public: // virtual methods

  /// Solve the linear system Ax = b
  virtual void apply ( std::shared_ptr<Matrix> A, Block * block) throw();

  /// Type of this solver
  virtual std::string type() const { return "jacobi"; }

#ifdef NEW_REFRESH  
  bool is_finest(Block * block) {return is_finest_(block); }
#endif  
  
protected: // virtual methods
  
#ifdef NEW_REFRESH  
  /// Whether Block is active
  virtual bool is_active_(Block * block) const
  {
    if (solve_type_ == solve_level) {
      return true;
    } else {
      return Solver::is_active_(block);
    }
  }

  /// Whether solution is defined on this Block
  virtual bool is_finest_(Block * block) const
  {
    if (solve_type_ == solve_level) {
      return true;
    } else {
      return Solver::is_finest_(block);
    }
  }
#endif
  
public: // methods

  /// Continue after refresh to perform Jacobi update
  void compute (Block * block);
  
protected: // methods

  /// Implementation of solver() for given precision 
  void apply_(Block * block);

  /// Refresh after computing
  void do_refresh_(Block * block);

  /// Allocate temporary Fields
  void allocate_temporary_(Field field, Block * block = NULL)
  {
    field.allocate_temporary(id_);
    field.allocate_temporary(ir_);
  }

  /// Dellocate temporary Fields
  void deallocate_temporary_(Field field, Block * block = NULL)
  {
    field.deallocate_temporary(id_);
    field.deallocate_temporary(ir_);
  }

  /// Return a pointer to the iteration counter on the block
  int * piter_(Block * block) {
    ScalarData<int> * scalar_data  = block->data()->scalar_data_int();
    ScalarDescr *     scalar_descr = cello::scalar_descr_int();
    return scalar_data->value(scalar_descr,i_iter_);
  }

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Matrix A for smoothing A*X = B
  std::shared_ptr<Matrix> A_;
  
  /// Field index for residual R
  int ir_;

  /// Field index for matrix diagonal D
  int id_;

  /// Weighting
  double w_;

  /// Scalar index for current iteration on a Block
  int i_iter_;
  
  /// Number of iterations
  int n_;

#ifdef NEW_REFRESH
  // Refresh after each smoothing
  int ir_smooth_;
#endif  
};

#endif /* ENZO_ENZO_SOLVER_JACOBI_HPP */

