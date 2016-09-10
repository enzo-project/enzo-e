// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialPm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

#ifndef ENZO_ENZO_INITIAL_PM_HPP
#define ENZO_ENZO_INITIAL_PM_HPP

class Mask;
class Parameters;

class EnzoInitialPm : public Initial {

  /// @class    EnzoInitialPm
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

public: // interface

  EnzoInitialPm (Parameters * parameters,
		 const std::string parameter_name,
		 int               init_cycle,
		 double            init_time,
		 std::string       field,
		 double            mpp,
		 int               level) throw ()
    : Initial(init_cycle, init_time),
      field_(field),
      mpp_(mpp),
      level_(level),
      mask_()
  {
    if (parameters) {
      mask_ = Mask::create (parameters->param(parameter_name),parameters);
    } else {
      mask_ = NULL;
    } 
  }
  
  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialPm);

  /// CHARM++ migration constructor
  EnzoInitialPm(CkMigrateMessage *m)
    : Initial (m),
      mpp_(0.0),
      level_(0),
      mask_(NULL)
  {  }

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

protected: // functions

  /// Initial particle positions are a uniform array if mpp_ <= 0
  void uniform_placement_ (Block * block, Field field, Particle particle);

  /// Initial particle positions are random based on local density
  void density_placement_ (Block * block, Field field, Particle particle);

private: // attributes

  /// Field to use for initial particle placement--default "density"
  std::string field_;

  /// Mass per particle--place on average one particle for each mpp grams
  double mpp_;

  /// Level corresponding to placing one particle per cell
  int level_;

  /// To define cloud extents
  Mask * mask_;

};

#endif /* ENZO_ENZO_INITIAL_PM_HPP */

