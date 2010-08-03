// $Id: user_UserTimestep.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     user_UserTimestep.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    Declaration for the UserTimestep component

#ifndef USER_USER_TIMESTEP_HPP
#define USER_USER_TIMESTEP_HPP


class UserTimestep {

  /// @class    UserTimestep
  /// @ingroup  User
  /// @brief    Encapsulate determination of timestep

public: // interface

  /// Create a new UserTimestep
  UserTimestep() throw()
  {};

public: // virtual functions

  /// Perform any method-specific initialization

  virtual void initialize (DataDescr * data_descr) throw() {} ;

  /// Perform any timestep-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  virtual void finalize (DataDescr * data_descr) throw(){};

  /// Initialize variables that may change for each block.  Called once per
  /// block per timestep.

  virtual void initialize_block (DataBlock * data_block) throw(){};

  /// Finalize after a timestep, e.g. to deallocate any
  /// dynamically-allocated variables

  virtual void finalize_block (DataBlock * data_block) throw(){};

  /// Compute the timestep for the block

  virtual double compute_block( DataBlock * data_block ) throw() = 0; 

};

#endif /* USER_USER_TIMESTEP_HPP */
