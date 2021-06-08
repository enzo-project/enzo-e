// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MaskExpr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the MaskExpr class
///

#ifndef PROBLEM_MASK_EXPR_HPP
#define PROBLEM_MASK_EXPR_HPP

#include "_parameters.hpp" // require full declaration of Expression class

class MaskExpr : public Mask {

  /// @class    MaskExpr
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Default Constructor
  MaskExpr() throw();

  /// Main Constructor
  MaskExpr(Parameters * parameters,
	   const std::string &parameter_name,
	   int parameter_index = -1) throw();

  /// Destructor
  virtual ~MaskExpr() throw() 
  { };

  /// Copy constructor
  MaskExpr(const MaskExpr & mask) = default;

  /// Assignment operator
  MaskExpr & operator= (const MaskExpr & mask) = default;

  /// Clone the object
  virtual std::shared_ptr<Mask> make_clone() const
  { return std::make_shared<MaskExpr> (*this); }

  PUPable_decl(MaskExpr);

  MaskExpr(CkMigrateMessage *m)
    : MaskExpr()
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Evaluate mask at a point
  virtual bool evaluate (double t, double x, double y, double z) const;

  /// Return mask values in an array
  virtual void evaluate (bool * mask, double t,
			 int ndx, int nx, double * x,
			 int ndy, int ny, double * y,
			 int ndz, int nz, double * z) const;

private: // attributes

  // NOTE: change pup() function whenever attributes change

  Expression expr_;

};

#endif /* PROBLEM_MASK_EXPR_HPP */

