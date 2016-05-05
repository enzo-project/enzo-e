// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialPm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

#ifndef ENZO_ENZO_INITIAL_PM_HPP
#define ENZO_ENZO_INITIAL_PM_HPP

class EnzoInitialPm : public Initial {

  /// @class    EnzoInitialPm
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

public: // interface

  EnzoInitialPm (int init_cycle, double init_time,
		 std::string field, double mpp) throw ()
    : Initial(init_cycle, init_time),
      field_(field),
      mpp_(mpp)
  { 
  }

  
  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialPm);

  /// CHARM++ migration constructor
  EnzoInitialPm(CkMigrateMessage *m) : Initial (m) {}

  /// Destructor
  virtual ~EnzoInitialPm() throw()
  {} ;

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  (
   Block * block,
   const FieldDescr * field_descr,
   const ParticleDescr * particle_descr,
   const Hierarchy * hierarchy
   ) throw();

private: // attributes

  /// Field to use for initial particle placement--default "density"
  std::string field_;

  /// Mass per particle--place on average one particle for each mpp grams
  double mpp_;

};

#endif /* ENZO_ENZO_INITIAL_PM_HPP */

