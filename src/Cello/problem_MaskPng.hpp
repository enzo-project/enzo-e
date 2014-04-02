// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MaskPng.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the MaskPng class
///

#ifndef PROBLEM_MASK_PNG_HPP
#define PROBLEM_MASK_PNG_HPP

class MaskPng : public Mask {

  /// @class    MaskPng
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  MaskPng() throw() 
  : Mask(),mask_(0), nx_(0),ny_(0), xm_(0),xp_(0),ym_(0),yp_(0)
  { };

  /// Destructor
  virtual ~MaskPng() throw() 
  { };

  /// Copy constructor
  MaskPng(const MaskPng & mask) throw() 
  {copy_(mask); }

  /// Assignment operator
  MaskPng & operator= (const MaskPng & mask) throw()
  {   copy_(mask); return *this; }


  /// Clone the object
  virtual Mask * clone() const
  { return (new MaskPng(*this)); }

  MaskPng(std::string file_name, 
	  double xm, double xp,
	  double ym, double yp) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    Mask::pup(p);
    p | nx_;
    p | ny_;
    PUParray(p,mask_,nx_*ny_);
    p | xm_;
    p | xp_;
    p | ym_;
    p | yp_;
    WARNING("MaskPng::pup()","UNFINISHED");
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

  void copy_(const MaskPng & mask) throw();

private: // attributes

  bool * mask_;
  /// Size of the image
  int nx_;
  int ny_;

  double xm_, xp_, ym_, yp_;

};

#endif /* PROBLEM_MASK_PNG_HPP */

