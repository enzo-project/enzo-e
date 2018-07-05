// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPmDeposit.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2016-04-25
/// @brief    [\ref Enzo] Implementation of mass-deposition for PM particle-mesh 

#ifndef ENZO_ENZO_METHOD_PM_DEPOSIT_HPP
#define ENZO_ENZO_METHOD_PM_DEPOSIT_HPP

class EnzoMethodPmDeposit : public Method {

  /// @class    EnzoMethodPmDeposit
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Declare Enzo's Particle-mesh method class

public: // interface

  /// Create a new EnzoMethodPmDeposit object
  EnzoMethodPmDeposit(double alpha = 0.0);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodPmDeposit);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodPmDeposit (CkMigrateMessage *m)
    : Method (m),
      alpha_(0.0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "pm_deposit"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

protected: // methods

    /// Allocate temporary Fields
  void allocate_temporary_(Block * block)
  {
    if (block->is_leaf()) {
      Field field = block->data()->field();
    }
  }

  /// Dellocate temporary Fields
  void deallocate_temporary_(Block * block)
  {
    if (block->is_leaf()) {
      Field field = block->data()->field();
    }
  }

protected: // attributes 

  /// Time parameter for when to deposit mass: time + alpha*dt
  double alpha_;
};

#endif /* ENZO_ENZO_METHOD_PM_DEPOSIT_HPP */
