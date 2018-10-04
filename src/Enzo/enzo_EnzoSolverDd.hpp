// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverDd.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2018-10-01
/// @brief    [\ref Enzo] Declaration of EnzoSolverDd
///
/// Domain decomposition solver

#ifndef ENZO_ENZO_SOLVER_DD_HPP
#define ENZO_ENZO_SOLVER_DD_HPP

class EnzoSolverDd : public Solver {

  /// @class    EnzoSolverDd
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Multigrid on the root-level grid using Mg0, then
  /// BiCgStab in overlapping subdomains defined by root-level Blocks.
  /// An optional final Jacobi step can be applied to smooth the solution
  /// along subdomain boundaries.

public: // interface

  /// Create a new EnzoSolverDd object
  EnzoSolverDd
  (std::string name,
   std::string field_x,
   std::string field_b,
   int index_solve_coarse,
   int index_solve_domain,
   int index_solve_smooth,
   int monitor_iter,
   int restart_cycle);

  EnzoSolverDd() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverDd);
  
  /// Charm++ PUP::able migration constructor
  EnzoSolverDd (CkMigrateMessage *m)
    :  Solver(m),
       A_(NULL),
       index_solve_coarse_(-1),
       index_solve_domain_(-1),
       index_solve_smooth_(-1),
       mx_(0),my_(0),mz_(0),
       gx_(0),gy_(0),gz_(0)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    // NOTE: change this function whenever attributes change

    TRACEPUP;

    Solver::pup(p);

    //    p | A_;
    p | index_solve_coarse_;
    p | index_solve_domain_;
    p | index_solve_smooth_;

    p | mx_;
    p | my_;
    p | mz_;
    p | gx_;
    p | gy_;
    p | gz_;

  }

  /// Solve the linear system 
  virtual void apply ( std::shared_ptr<Matrix> A, Block * block) throw();

  /// Type of this solver
  virtual std::string type() const { return "dd"; }

  /// Allocate temporary Fields
  void allocate_temporary_(Block * block)
  {
  }
	      
  /// Dellocate temporary Fields

  void deallocate_temporary_(Block * block)
  {
  }

  /// End of solver
  void end(Block* block) throw();
  
protected: // methods

protected: // attributes

  /// Matrix
  std::shared_ptr<Matrix> A_;

  /// Indices for coarse solver, domain solver, and smoother

  int index_solve_coarse_;
  int index_solve_domain_;
  int index_solve_smooth_;
  
  /// Block field attributes
  int mx_,my_,mz_;
  int gx_,gy_,gz_;

};

#endif /* ENZO_ENZO_SOLVER_GRAVITY_DD_HPP */
