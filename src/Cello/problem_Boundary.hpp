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
  Boundary() throw() 
  : axis_(axis_all), face_(face_all), mask_(nullptr)
  {
    for (int axis=0; axis<3; axis++) {
      periodicity_[axis] = false;
    }
  }

  /// Create a new Boundary
  Boundary(axis_enum axis, face_enum face, std::shared_ptr<Mask> mask) throw() 
  : axis_(axis), face_(face), mask_(mask)
  {
    for (int axis=0; axis<3; axis++) {
      periodicity_[axis] = false;
    }
  }

  /// Destructor
  virtual ~Boundary() throw()
  {
  }

  /// Charm++ PUP::able declarations
  PUPable_abstract(Boundary);

  /// CHARM++ migration constructor for PUP::able
  Boundary (CkMigrateMessage *m)
    : PUP::able(m),
       axis_(axis_all), face_(face_all), mask_(nullptr)
  {
    for (int axis=0; axis<3; axis++) {
      periodicity_[axis] = false;
    }
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; 
    PUP::able::pup(p); 
    int axis=axis_;
    p | axis;
    axis_ = (axis_enum)axis;
    int face = face_;
    p | face;
    face_ = (face_enum)face;
    static bool warn[CONFIG_NODE_SIZE] = {false};
    const int in = cello::index_static();
    if (! warn[in]) {
      WARNING ("Boundary::pup()", "Skipping p | mask_");
      warn[in] = true;
    }
    PUParray(p,periodicity_,3);
  };

public: // virtual functions

  /// Enforce boundary conditions

  virtual void enforce (Block   * block,
			face_enum face = face_all,
			axis_enum axis = axis_all) const throw() = 0;

  /// Return which faces are periodic
  void periodicity(bool p3[3]) const throw() {
    for (int axis=0; axis<3; axis++) {
      if (periodicity_[axis]) p3[axis] = true;
    }
  }

protected: // protected functions

  /// Return whether the boundary condition applies for the given axis and face
  bool applies_(axis_enum axis, face_enum face) const throw ()
  { return ((axis_==axis_all || axis==axis_) && 
	    (face_==face_all || face==face_)); }

protected: // protected attributes

  /// Axis mask for boundary conditions (may be axis_all)
  axis_enum axis_; 

  /// Face mask for boundary conditions (may be face_all)
  face_enum face_;

  /// Mask object for boundary conditions (NULL == true)
  std::shared_ptr<Mask> mask_;

  /// Periodicity of boundary condition faces
  bool periodicity_[3];

};

#endif /* PROBLEM_BOUNDARY_HPP */
