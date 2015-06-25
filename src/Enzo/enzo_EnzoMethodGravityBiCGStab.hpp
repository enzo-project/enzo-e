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
  EnzoMethodGravityBiCGStab(const FieldDescr * field_descr,
			    int iter_max, double res_tol);

  EnzoMethodGravityBiCGStab() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityBiCGStab);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityBiCGStab (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Solve for the gravitational potential
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "gravity_bicgstab"; }

protected: // methods

  template <class T>
  void compute_ (Block * block,
		 T * density,   int md3[3], int nd3[3],
		 T * potential, int mp3[3], int np3[3]) const throw();

protected: // attributes

  /// Maximum number of BiCGStab iterations
  int iter_max_;

  /// Convergence tolerance on the residual 
  double res_tol_;

};

#endif /* ENZO_ENZO_METHOD_GRAVITY_BICGSTAB_HPP */
