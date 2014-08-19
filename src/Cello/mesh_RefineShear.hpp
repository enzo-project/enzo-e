// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineShear.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 21 16:23:57 PDT 2014
/// @brief    [\ref Mesh] Declaration of the RefineShear class
///

#ifndef MESH_REFINE_SHEAR_HPP
#define MESH_REFINE_SHEAR_HPP

class RefineShear : public Refine {

  /// @class    RefineShear
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  RefineShear(const FieldDescr * field_descr,
	      double shear_min_refine,
	      double shear_max_coarsen,
	      std::string output) throw();

  /// default constructor
  // RefineShear () throw() : Refine() {};

  PUPable_decl(RefineShear);

  RefineShear(CkMigrateMessage *m) : Refine (m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    Refine::pup(p);
  }

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (CommBlock * comm_block,
		     const FieldDescr * field_descr) throw();

  virtual std::string name () const { return "shear"; };

private: // functions

  template <class T>
  void evaluate_block_(const T * u,
		       const T * v,
		       const T * w,
		       T * output,
		       int ndx, int ndy, int ndz,
		       int nx, int ny, int nz,
		       int gx, int gy, int gz,
		       bool *any_refine,
		       bool * all_coarsen, 
		       int rank);

};

#endif /* MESH_REFINE_SHEAR_HPP */

