// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Mask.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    [\ref Problem] Declaration of the Mask class
///

#ifndef PROBLEM_MASK_HPP
#define PROBLEM_MASK_HPP

class Param;
class Parameters;

class Mask : public PUP::able {

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
  virtual std::shared_ptr<Mask> make_clone() const = 0;

  /// Create a new object of the appropriate subclass
  static std::shared_ptr<Mask> create(Parameters * parameters,
				      const std::string &parameter_name,
				      int parameter_index = -1);

  /// CHARM++ PUP::able declaration
  PUPable_abstract(Mask);

  /// CHARM++ migration constructor for PUP::able
  Mask (CkMigrateMessage *m) : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    PUP::able::pup(p); 
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

