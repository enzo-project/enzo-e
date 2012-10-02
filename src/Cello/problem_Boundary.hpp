// See LICENSE_CELLO file for license and copyright information

/// @file     method_Boundary.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Method] Declaration for the Boundary component

#ifndef METHOD_BOUNDARY_HPP
#define METHOD_BOUNDARY_HPP


#ifdef CONFIG_USE_CHARM
class Boundary : public PUP::able 
#else
class Boundary 
#endif
{

  /// @class    Boundary
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate a boundary conditions generator

public: // interface

  /// Create a new Boundary
  Boundary() throw() {}

  /// Destructor
  virtual ~Boundary() throw() {}

#ifdef CONFIG_USE_CHARM

  /// Charm++ PUP::able declarations
  PUPable_abstract(Boundary);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#endif


public: // virtual functions

  /// Enforce boundary conditions

  virtual void enforce (const FieldDescr * field_descr,
			Block * block,
			face_enum face = face_all,
			axis_enum axis = axis_all) const throw() = 0;

  /// Whether boundary conditions are periodic (handled by ghost refresh)
  virtual bool is_periodic() const throw() = 0;

};

#endif /* METHOD_BOUNDARY_HPP */
