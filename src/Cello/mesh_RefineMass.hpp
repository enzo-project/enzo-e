// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineMass.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the RefineMass class
///

#ifndef MESH_REFINE_MASS_HPP
#define MESH_REFINE_MASS_HPP

class RefineMass : public Refine {

  /// @class    RefineMass
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  RefineMass(double min_mass,
	     double level_exponent,
	     double min_overdensity,
	     double root_cell_volume) throw();

#ifdef CONFIG_USE_CHARM

  /// default constructor
  RefineMass () throw() : Refine() {};

  PUPable_decl(RefineMass);

  RefineMass(CkMigrateMessage *m) : Refine (m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    Refine::pup(p);
    p | min_;
    p | level_exponent_;
    p | min_overdensity_;
  }
#endif

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (FieldBlock * field_block,
		     const FieldDescr * field_descr) throw();

private:

  /// Minimum allowed mass before refinement kicks in
  double min_;

  /// Minimum allowed mass before refinement kicks in
  double level_exponent_;

  /// Minimum allowed mass overdensity before refinement kicks in
  double min_overdensity_;
};

#endif /* MESH_REFINE_MASS_HPP */

