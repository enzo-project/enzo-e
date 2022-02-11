// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Ordering.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-02-01
/// @brief    [\ref Mesh] Declaration of the Ordering class

#ifndef MESH_ORDERING_HPP
#define MESH_ORDERING_HPP

class Ordering : public PUP::able {

  /// @class    Ordering
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Ordering() throw()
  { }

  /// CHARM++ PUP::able declaration
  PUPable_abstract(Ordering);


  /// CHARM++ migration constructor for PUP::able
  Ordering (CkMigrateMessage *m) :
    PUP::able(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  {
    TRACEPUP;
    PUP::able::pup(p); 
  }

public:  // virtual methods

  /// Return the index of the given Block
  virtual int index (Block * block) const = 0;

  /// Return the number of indices in the Ordering
  virtual int num_indices () const = 0;

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* MESH_ORDERING_HPP */

