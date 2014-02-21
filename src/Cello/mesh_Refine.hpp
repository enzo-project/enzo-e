// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Refine.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the Refine Base class
///

#ifndef MESH_REFINE_HPP
#define MESH_REFINE_HPP

class Refine : public PUP::able 
{

  /// @class    Refine
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  Refine() throw() {};

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
  
  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (CommBlock * comm_block,
		     const FieldDescr * field_descr) throw ()
  {
    ERROR("Refine::apply",
	 "'Abstract' virtual function Refine::apply() should not be called");
    return 0;
  };

  virtual std::string name () const { return "unknown"; }

};

#endif /* MESH_REFINE_HPP */

