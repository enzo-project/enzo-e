// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Refine.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the Refine Base class
///

#ifndef MESH_REFINE_HPP
#define MESH_REFINE_HPP

class Refine 
#ifdef CONFIG_USE_CHARM
  : public PUP::able 
#endif
{

  /// @class    Refine
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Refine() throw() {};

#ifdef CONFIG_USE_CHARM
  /// CHARM++ PUP::able declaration
  PUPable_decl(Refine);

  /// CHARM++ migration constructor for PUP::able

  Refine (CkMigrateMessage *m) : PUP::able(m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    PUP::able::pup(p);
    // NOTE: change this function whenever attributes change
  }
#endif
  
  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (FieldBlock * field_block,
		     const FieldDescr * field_descr) throw ()
  {
    ERROR("Refine::apply",
	 "'Abstract' virtual function Refine::apply() should not be called");
  };

};

#endif /* MESH_REFINE_HPP */

