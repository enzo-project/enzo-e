// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Mask.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the Mask class
///

#ifndef PROBLEM_MASK_HPP
#define PROBLEM_MASK_HPP

class Mask {

  /// @class    Mask
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  Mask() throw() {};

  /// Destructor
  virtual ~Mask() throw() {};

  /// Copy constructor
  Mask(const Mask & Mask) throw() {};

  /// Assignment operator
  Mask & operator= (const Mask & Mask) throw() {return *this; }

  /// Clone the object
  virtual Mask * clone() const = 0;

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
  }

  /// Evaluate mask at a point
  virtual bool evaluate (double t, double x, double y, double z) const = 0;

  /// Return mask values in an array
  virtual void evaluate (bool * mask, double t,
			 int ndx, int nx, double * x,
			 int ndy, int ny, double * y,
			 int ndz, int nz, double * z) const = 0;

  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* PROBLEM_MASK_HPP */

