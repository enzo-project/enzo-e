// See LICENSE_CELLO file for license and copyright information

/// @file     problem_BoundaryValue.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-04-01
/// @brief    [\ref Problem] Declaration for the BoundaryValue component

#ifndef PROBLEM_BOUNDARY_VALUE_HPP
#define PROBLEM_BOUNDARY_VALUE_HPP

class BoundaryValue : public Boundary
{

  /// @class    BoundaryValue
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Encapsulate a BoundaryValue conditions generator

public: // interface

  /// Create a new BoundaryValue
  BoundaryValue() throw() 
  : Boundary (), value_(0) 
  {  }

  /// Create a new BoundaryValue
  BoundaryValue(axis_enum axis, face_enum face, Value * value, 
		std::vector<std::string> field_list) throw() 
    : Boundary(axis,face,0), value_(value), field_list_(field_list)
  { }

  /// Destructor
  virtual ~BoundaryValue() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(BoundaryValue);

  BoundaryValue(CkMigrateMessage *m) : Boundary (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  {
    // NOTE: change this function whenever attributes change
    Boundary::pup(p); 
    TRACEPUP; 

    p | *value_;
    p | field_list_;
  };

public: // virtual functions

  /// Enforce BoundaryValue conditions

  virtual void enforce (Block   * block,
			face_enum face = face_all,
			axis_enum axis = axis_all) const throw();

protected: // functions

  template <class T>
  void copy_(T * field, double * value,
	     int ndx, int ndy, int ndz,
	     int nx,  int ny,  int nz,
	     int ix0, int iy0, int iz0) const throw ();

protected: // attributes

  Value * value_;
  std::vector<std::string> field_list_;

};

#endif /* PROBLEM_BOUNDARY_VALUE_HPP */
