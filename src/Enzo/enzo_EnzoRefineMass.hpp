// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRefineMass.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-08-14
/// @brief    [\ref Enzo] Declaration of the EnzoRefineMass class
///

#ifndef ENZO_REFINE_MASS_HPP
#define ENZO_REFINE_MASS_HPP

class EnzoRefineMass : public Refine {

  /// @class    EnzoRefineMass
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] 

public: // interface

  /// Constructor
  EnzoRefineMass(double min_refine,
		 double max_coarsen,
		 int    max_level,
		 bool include_ghosts,
		 std::string output,
		 std::string name,
		 std::string mass_type,
		 double level_exponent) throw();

  /// default constructor
  // EnzoRefineMass () throw() : Refine() {};

  PUPable_decl(EnzoRefineMass);

  EnzoRefineMass(CkMigrateMessage *m)
    : Refine (m),
      name_(""),
      mass_ratio_(0.0),
      level_exponent_(0.0)
  { }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    Refine::pup(p);
    p | name_;
    p | mass_ratio_;
    p | level_exponent_;
  }

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (Block * block) throw();

  virtual std::string name () const { return "mass"; };

private:

  /// Field containing density to compare against
  std::string name_;

  /// Ratio of type of mass (baryon or cdm) to total if cosmology, else 0.0
  double mass_ratio_;

  /// Level expontent
  double level_exponent_;
};

#endif /* ENZO_REFINE_MASS_HPP */

