// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineDensity.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-07-31
/// @brief    [\ref Mesh] Declaration of the RefineDensity class
///

#ifndef MESH_REFINE_DENSITY_HPP
#define MESH_REFINE_DENSITY_HPP

class RefineDensity : public Refine {

  /// @class    RefineDensity
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 
  ///
  /// RefineDensity refines a block if any density field values are
  /// greater than a maximum density threshold, or coarsens if all
  /// density field values are less than a minimum density threshold

public: // interface

  /// Constructor
  RefineDensity(double min_refine,
		double max_coarsen,
		int    max_level,
		bool include_ghosts,
		std::string output) throw();

  PUPable_decl(RefineDensity);

  RefineDensity(CkMigrateMessage *m) : Refine (m) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    Refine::pup(p);
  }

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (Block * block) throw();

  virtual std::string name () const { return "density"; };

private: // functions

  template <class T>
  int apply_ (const T * array,
	      int mx, int my, int mz,
	      int gx, int gy, int gz) const throw ();

};

#endif /* MESH_REFINE_DENSITY_HPP */

