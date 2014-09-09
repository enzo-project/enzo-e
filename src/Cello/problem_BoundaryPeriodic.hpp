// See LICENSE_CELLO file for license and copyright information

/// @file     problem_BoundaryPeriodic.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-05-06
/// @brief    [\ref Problem] Declaration for the BoundaryPeriodic class

#ifndef PROBLEM_BOUNDARY_PERIODIC_HPP
#define PROBLEM_BOUNDARY_PERIODIC_HPP

class BoundaryPeriodic : public Boundary
{

  /// @class    BoundaryPeriodic
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Periodic boundary conditions

public: // interface

  /// Create a new BoundaryPeriodic
  BoundaryPeriodic() throw() 
  {
    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
	periodicity_[axis][face] = true;
      }
    }
  }

  /// Destructor
  virtual ~BoundaryPeriodic() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(BoundaryPeriodic);

  BoundaryPeriodic(CkMigrateMessage *m) : Boundary (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; 
    PUP::able::pup(p); 
  };

public: // virtual functions

  /// Enforce BoundaryPeriodic conditions

  virtual void enforce (CommBlock * block,
			face_enum face = face_all,
			axis_enum axis = axis_all) const throw()
  { };

};

#endif /* PROBLEM_BOUNDARY_PERIODIC_HPP */
