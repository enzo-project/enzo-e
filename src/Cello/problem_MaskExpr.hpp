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
  MaskExpr() throw() 
  : Mask() , param_(0)
  { };

  /// Destructor
  virtual ~MaskExpr() throw() 
  { };

  /// Copy constructor
  MaskExpr(const MaskExpr & mask) throw() 
  {copy_(mask); }

  /// Assignment operator
  MaskExpr & operator= (const MaskExpr & mask) throw()
  {   copy_(mask); return *this; }


  /// Clone the object
  virtual Mask * clone() const
  { return (new MaskExpr(*this)); }

  MaskExpr(Param * param) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    Mask::pup(p);

    WARNING("MaskExpr::pup()","UNFINISHED");
    // NOTE: change this function whenever attributes change
  }

  /// Evaluate mask at a point
  virtual bool evaluate (double t, double x, double y, double z) const;

  /// Return mask values in an array
  virtual void evaluate (bool * mask, double t,
			 int ndx, int nx, double * x,
			 int ndy, int ny, double * y,
			 int ndz, int nz, double * z) const;
  
private: // functions

  void copy_(const MaskExpr & mask) throw();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  Param * param_;


};

#endif /* PROBLEM_MASK_EXPR_HPP */

