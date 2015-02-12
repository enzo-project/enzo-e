// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineMask.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the RefineMask class
///

#ifndef MESH_REFINE_MASK_HPP
#define MESH_REFINE_MASK_HPP

class Value;
class Parameters;

class RefineMask : public Refine {

  /// @class    RefineMask
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  RefineMask(Parameters * parameters,
	     const std::string parameter_name,
	     std::string output) throw();

  /// default constructor
  // RefineMask () throw() : Refine() {};

  PUPable_decl(RefineMask);

  RefineMask(CkMigrateMessage *m) : Refine (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (Block            * block,
		     const FieldDescr * field_descr) throw();

  virtual std::string name () const { return "mask"; };

private: // functions

private: // attributes

  Value * value_;
};

#endif /* MESH_REFINE_MASK_HPP */

