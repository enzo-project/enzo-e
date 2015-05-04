// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeSmoothJacobi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-04-30 18:45:58
/// @brief    [\ref Enzo] Declaration of the EnzoComputeSmoothJacobi class

#ifndef ENZO_ENZO_COMPUTE_SMOOTH_JACOBI_HPP
#define ENZO_ENZO_COMPUTE_SMOOTH_JACOBI_HPP

class EnzoComputeSmoothJacobi : public Compute {

  /// @class    EnzoComputeSmoothJacobi
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoComputeSmoothJacobi(std::string x_field,
			  std::string r_field,
			  std::string d_field,
			  double weight,
			  FieldDescr * field_descr) throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoComputeSmoothJacobi);

  /// Charm++ PUP::able migration constructor
  EnzoComputeSmoothJacobi (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { 
    TRACEPUP;
    Compute::pup(p);

    p | A_;
    p | i_x_;
    p | i_r_;
    p | i_d_;
  }
  
public: // virtual functions

  virtual void compute ( Block * block) throw(); 

private: // functions

  /// Implementation of compute() for given precision 
  template <typename T>
  void compute_(Block * block);

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Matrix A for smoothing A*X = B
  Matrix * A_;
  
  /// Field index for vector X
  int i_x_;

  /// Field index for residual R
  int i_r_;

  /// Field index for matrix diagonal D
  int i_d_;

  /// Weighting
  double w_;
};

#endif /* ENZO_ENZO_COMPUTE_SMOOTH_JACOBI_HPP */

