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
  Refine(double min_refine,
	 double max_coarsen,
	 std::string output) throw()
    : min_refine_(min_refine),
      max_coarsen_(max_coarsen),
      output_(output)
  {};

  /// CHARM++ PUP::able declaration
  PUPable_decl(Refine);

  /// CHARM++ migration constructor for PUP::able

  Refine (CkMigrateMessage *m) : PUP::able(m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    PUP::able::pup(p);
    // NOTE: change this function whenever attributes change
    p | min_refine_;
    p | max_coarsen_;
    p | output_;
  }
  
  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (Block * block,
		     const FieldDescr * field_descr) throw ()
  {
    ERROR("Refine::apply",
	 "'Abstract' virtual function Refine::apply() should not be called");
    return 0;
  };

  virtual std::string name () const { return "unknown"; }

  void * initialize_output_(FieldBlock * field_block);

protected:

  /// Minimum allowed value before refinement kicks in
  double min_refine_;

  /// Maximum allowed value before coarsening is allowed
  double max_coarsen_;

  /// Field name to write refinement field to (-1 coarsen 0 same +1 refine)
  std::string output_;
  
};

#endif /* MESH_REFINE_HPP */

