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
  RefineMass(double min_refine,
	     double max_coarsen,
	     double level_exponent,
	     double root_cell_volume) throw();

  /// default constructor
  RefineMass () throw() : Refine() {};

  PUPable_decl(RefineMass);

  RefineMass(CkMigrateMessage *m) : Refine (m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    Refine::pup(p);
    p | min_refine_;
    p | max_coarsen_;
    p | level_exponent_;
  }

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (CommBlock * comm_block,
		     const FieldDescr * field_descr) throw();

  virtual std::string name () const { return "mass"; };

private:

  /// Minimum allowed mass before refinement kicks in
  double min_refine_;

  /// Minimum allowed mass before refinement kicks in
  double max_coarsen_;

  double level_exponent_;

};

#endif /* MESH_REFINE_MASS_HPP */

