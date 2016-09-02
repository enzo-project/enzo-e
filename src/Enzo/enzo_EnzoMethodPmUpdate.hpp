// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPmUpdate.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2016-04-25
/// @brief    [\ref Enzo] Implementation of particle updating for for PM method

#ifndef ENZO_ENZO_METHOD_PM_UPDATE_HPP
#define ENZO_ENZO_METHOD_PM_UPDATE_HPP

class EnzoMethodPmUpdate : public Method {

  /// @class    EnzoMethodPmUpdate
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] PM method particle update

public: // interface

  /// Create a new EnzoMethodPmUpdate object
  EnzoMethodPmUpdate(const FieldDescr * , 
		     const ParticleDescr *, double max_dt);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodPmUpdate);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodPmUpdate (CkMigrateMessage *m)
    : max_dt_(0.0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "pm_update"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // attributes

  double max_dt_;

};

#endif /* ENZO_ENZO_METHOD_PM_UPDATE_HPP */
