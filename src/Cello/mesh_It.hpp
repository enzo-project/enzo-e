// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_It.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Aug 25 15:03:38 PDT 2011
/// @brief    [\ref Mesh] Declaration of the It class

#ifndef MESH_IT_HPP
#define MESH_IT_HPP

template<class T>
class It {

  /// @class    It
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator abstract base class

public: // interface

  /// It constructor
  It () throw () : index1_(0) {};

  /// It destructor
  virtual ~It() throw()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | index1_;
  }

  /// Iterate through entities
  virtual T * operator++ () throw() = 0;

  /// Return whether the iteration is complete
  virtual bool done() const throw() = 0;


protected: // attributes

  /// Index of the current entity plus 1, or 0 if between iterations
  /// Always in the range 0 <= index1_ <= number of local patchs
  size_t index1_;
};

#endif /* MESH_IT_HPP */
