// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MaskExpr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the MaskExpr class
///

#ifndef PROBLEM_MASK_EXPR_HPP
#define PROBLEM_MASK_EXPR_HPP

class MaskExpr : public Mask {

  /// @class    MaskExpr
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  MaskExpr() throw() : Mask() { };

  /// Destructor
  ~MaskExpr() throw();

  /// Copy constructor
  MaskExpr(const MaskExpr & mask) throw();

  /// Assignment operator
  MaskExpr & operator= (const MaskExpr & mask) throw();

  MaskExpr(Parameters * parameters,
	   const std::string parameter_name,
	   int index_parameter,
	   double time,
	   int nx, double xmin, double xmax,
	   int ny=1, double ymin=0, double ymax=0,
	   int nz=1, double zmin=0, double zmax=0) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    Mask::pup(p);
    
    // NOTE: change this function whenever attributes change
  }

  virtual bool evaluate (int ix, int iy, int iz) const;
  virtual void evaluate (bool * mask, int ndx, int ndy, int ndz) const;
  
private: // functions

  void copy_(const MaskExpr & mask) throw();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  Parameters * parameters_;
  std::string parameter_name_;
  int index_parameter_;
  int nx_,ny_,nz_;
  double time_;
  double xmin_,ymin_,zmin_;
  double xmax_,ymax_,zmax_;


};

#endif /* PROBLEM_MASK_EXPR_HPP */

