// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Boundary.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Problem] Declaration for the Boundary component

#ifndef PROBLEM_BOUNDARY_HPP
#define PROBLEM_BOUNDARY_HPP


class Boundary : public PUP::able 
{

  /// @class    Boundary
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Encapsulate a boundary conditions generator

public: // interface

  /// Create a new Boundary
  Boundary() throw() {}

  /// Destructor
  virtual ~Boundary() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_abstract(Boundary);

  /// CHARM++ migration constructor for PUP::able
  Boundary (CkMigrateMessage *m) : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; PUP::able::pup(p); };

public: // virtual functions

  /// Enforce boundary conditions

  virtual void enforce (const FieldDescr * field_descr,
			CommBlock * block,
			face_enum face = face_all,
			axis_enum axis = axis_all) const throw() = 0;

};

#endif /* PROBLEM_BOUNDARY_HPP */
