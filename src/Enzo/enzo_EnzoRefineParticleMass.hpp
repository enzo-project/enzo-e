// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRefineParticleMass.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-05-16
/// @brief    [\ref Mesh] Declaration of the EnzoRefineParticleMass class
///

#ifndef ENZO_REFINE_PARTICLE_MASS_HPP
#define ENZO_REFINE_PARTICLE_MASS_HPP

class EnzoRefineParticleMass : public Refine {

  /// @class    EnzoRefineParticleMass
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  EnzoRefineParticleMass
  (double min_refine,
   double max_coarsen,
   int    max_level,
   bool   include_ghosts,
   std::string store,
   double level_exponent
   ) throw();

  PUPable_decl(EnzoRefineParticleMass);

  EnzoRefineParticleMass(CkMigrateMessage *m) :
    Refine (m),
    level_exponent_(0.0) {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    Refine::pup(p);
    p | level_exponent_;
  }

  /// Evaluate the refinement criteria, updating the refinement field
  virtual int apply (Block * block) throw();

  virtual std::string name () const { return "particle_mass"; };

private:

  double level_exponent_;

};

#endif /* ENZO_REFINE_PARTICLE_MASS_HPP */

