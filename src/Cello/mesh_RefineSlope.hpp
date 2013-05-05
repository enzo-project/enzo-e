// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineSlope.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Mesh] Declaration of the RefineSlope class
///

#ifndef MESH_REFINE_SLOPE_HPP
#define MESH_REFINE_SLOPE_HPP

class RefineSlope : public Refine {

  /// @class    RefineSlope
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  RefineSlope(const FieldDescr * field_descr,
	      double slope_min_refine,
	      double slope_max_coarsen,
	      std::vector<std::string> field_name_list) throw();

#ifdef CONFIG_USE_CHARM

  /// default constructor
  RefineSlope () throw() : Refine() {};

  PUPable_decl(RefineSlope);

  RefineSlope(CkMigrateMessage *m) : Refine (m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    Refine::pup(p);
    p | slope_min_refine_;
    p | slope_max_coarsen_;
    p | field_id_list_;
  }
#endif
  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (CommBlock * comm_block,
		     const FieldDescr * field_descr) throw();

  virtual std::string name () const { return "slope"; };

private: // functions

  template <class T>
  void evaluate_block_(T * array, 
		       int nx, int ny, int nz,
		       int gx, int gy, int gz,
		       bool *any_refine,
		       bool * all_coarsen, 
		       int rank, 
		       double * h3, const int * d3);

private: // attributes

  /// Minimum allowed slope before refinement kicks in
  double slope_min_refine_;

  /// Maximum allowed slope before coarsening is allowed
  double slope_max_coarsen_;

  /// List of field id's
  std::vector <int> field_id_list_;
  
  /// TEMPORARY
  bool debug_;
};

#endif /* MESH_REFINE_SLOPE_HPP */

