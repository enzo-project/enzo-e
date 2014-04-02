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
  BoundaryValue() throw() {}

  /// Destructor
  virtual ~BoundaryValue() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(BoundaryValue);

  BoundaryValue(CkMigrateMessage *m) : Boundary (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; PUP::able::pup(p); };

public: // virtual functions

  /// Enforce BoundaryValue conditions

  virtual void enforce (const FieldDescr * field_descr,
			CommBlock * block,
			face_enum face = face_all,
			axis_enum axis = axis_all) const throw() = 0;

};

#endif /* PROBLEM_BOUNDARY_VALUE_HPP */
