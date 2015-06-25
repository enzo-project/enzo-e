// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityBiCGStab.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     2014-10-21 17:25:40
/// @brief    [\ref Enzo] Declaration of EnzoMethodGravityBiCGStab
///
/// Bicongugate gradient stabilized method (BiCGStab) for solving 
/// for self-gravity on field data.

#ifndef ENZO_ENZO_METHOD_GRAVITY_BICGSTAB_HPP
#define ENZO_ENZO_METHOD_GRAVITY_BICGSTAB_HPP

class EnzoMethodGravityBiCGStab : public Method {

  /// @class    EnzoMethodGravityBiCGStab
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve self-gravity
  /// using the BiCGStab method.  Alone, this is more applicable to 
  /// smaller problems since the method doesn't scale as well as 
  /// some other methods (FFT, MG, etc.) for larger problems.  
  /// Alternately, a more scalable method may be combined as a 
  /// preconditioner for a robust and scalable overall method.

public: // interface

  /// Create a new EnzoMethodGravityBiCGStab object
  EnzoMethodGravityBiCGStab(const FieldDescr* field_descr,
			    int rank,
			    double grav_const,
			    int iter_max, 
			    double res_tol,
			    int monitor_iter,
			    bool is_singular,
			    bool diag_precon);

  EnzoMethodGravityBiCGStab() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityBiCGStab);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityBiCGStab(CkMigrateMessage* m) {}

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er& p) {

    // NOTE: change this function whenever attributes change
    TRACEPUP;

    Method::pup(p);

    p | A_;
    p | M_;
    p | is_singular_;
    p | rank_;
    p | grav_const_;
    p | iter_max_;
    p | res_tol_;
    p | monitor_iter_;
    p | rr0_;
    p | rr_min_;
    p | rr_max_;
    p | idensity_;
    p | ipotential_;
    p | ib_;
    p | ix_;
    p | ir_;
    p | id_;
    p | iy_;
    p | iz_;
    p | is_leaf_;

    p | nx_;
    p | ny_;
    p | nz_;

    p | mx_;
    p | my_;
    p | mz_;

    p | gx_;
    p | gy_;
    p | gz_;

    p | iter_;
  
    p | rz_;
    p | rz2_;
    p | dy_;
    p | bs_;
    p | bc_;
    p | id_refresh_matvec_;

  }

  
  /// Solve for the gravitational potential
  virtual void compute(Block* block) throw();

  virtual std::string name() throw () { return "gravity_bicgstab"; }

  /// Continuation after global reduction
  template<class T> void bicgstab_shift_1(EnzoBlock* enzo_block) throw();

  /// Continuation after global reduction
  template<class T> void bicgstab_loop_2(EnzoBlock* enzo_block) throw();

  /// Continuation after global reduction
  template<class T> void bicgstab_loop_4(EnzoBlock* enzo_block) throw();

  // /// Continuation after global reduction
  // template<class T> void bicgstab_shift_2(EnzoBlock* enzo_block) throw();

  /// Continuation after global reduction
  template<class T> void bicgstab_loop_6(EnzoBlock* enzo_block) throw();

  template<class T> void bicgstab_end(EnzoBlock* enzo_block, int retval) throw();

  /// Set rz_ by EnzoBlock after reduction
  void set_rz(long double rz) throw() { rz_ = rz; }

  /// Set rr_ by EnzoBlock after reduction
  void set_rr(long double rr) throw() { rr_ = rr; }

  /// Set rr_new_ by EnzoBlock after reduction
  void set_rz2(long double rz2) throw() { rz2_ = rz2; }

  /// Set dy_ by EnzoBlock after reduction
  void set_dy(long double dy) throw() { dy_ = dy; }

  /// Set bs_ (B sum) by EnzoBlock after reduction
  void set_bs(long double bs) throw() { bs_ = bs; }
  /// Set rs_ (R sum) by EnzoBlock after reduction
  void set_rs(long double rs) throw() { rs_ = rs; }
  /// Set xs_ (X sum) by EnzoBlock after reduction
  void set_xs(long double xs) throw() { xs_ = xs; }

  /// Set bc_ (B count) by EnzoBlock after reduction
  void set_bc(long double bc) throw() { bc_ = bc; }

  /// Set iter_ by EnzoBlock after reduction
  void set_iter(int iter) throw() { iter_ = iter; }

protected: // methods

  void monitor_output_(EnzoBlock * enzo_block) throw();

  template<class T> void compute_(EnzoBlock * enzo_block) throw();

  template<class T> void bicgstab_begin_1_() throw();

  void bicgstab_exit_() throw();

  /// Compute local contribution to inner-product X*Y
  template<class T> long double dot_(const T* X, const T* Y) const throw();

  template<class T> void zaxpy_(T* Z, double a, const T* X, const T* Y) const throw();
  
  /// Compute local sum of vector elements X_i
  template<class T> long double sum_(const T* X) const throw();

  /// scale the vector by the given scalar Y = a*X
  template<class T> void scale_(T* Y, T a, const T* X) const throw();

  /// return the number of elements of the vector X
  template<class T> int count_(T* X) const throw();
  
  /// Shift the vector X by a scalar multiple of Y
  /// NOTE includes ghost zones since performed after ghost refresh
  template<class T> void shift_(T* X, const T a, const T* Y) const throw();

  /// Set whether current Block is a leaf--if not don't touch data
  void set_leaf(Block* block) throw() { is_leaf_ = block->is_leaf(); }

protected: // attributes

  /// Matrix
  Matrix* A_;

  /// Preconditioner
  Matrix* M_;

  /// Whether you need to subtract of the nullspace of A from b, e.g. fully
  /// periodic or Neumann problems
  bool is_singular_;

  /// Dimensionality of the problem
  int rank_;

  /// Gas constant, e.g. 6.67384e-8 [cgs]
  double grav_const_;

  /// Maximum number of BiCGStab iterations
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

  /// BiCGStab vector id's
  int ib_;
  int ix_;
  int ir_;
  int id_;
  int iy_;
  int iz_;

  /// Whether current block is a leaf
  bool is_leaf_;

  /// Block field attributes
  int nx_, ny_, nz_;
  int mx_, my_, mz_;
  int gx_, gy_, gz_;

  /// Current BiCGStab iteration
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

#endif /* ENZO_ENZO_METHOD_GRAVITY_BICGSTAB_HPP */
