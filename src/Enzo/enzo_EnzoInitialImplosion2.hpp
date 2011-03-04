// $Id: enzo_EnzoInitialImplosion2.hpp 1877 2010-11-30 01:20:27Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_INITIAL_IMPLOSION2_HPP
#define ENZO_ENZO_INITIAL_IMPLOSION2_HPP

/// @file     enzo_EnzoInitialImplosion2.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

class EnzoInitialImplosion2 : public Initial {

  /// @class    EnzoInitialImplosion2
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

public: // interface

  /// Constructor
  EnzoInitialImplosion2(Monitor   * monitor,
                        EnzoDescr * enzo) throw();

  /// Initialize the block

  virtual void compute (DataBlock * data_block) throw();

private: // attributes

  EnzoDescr * enzo_;
  
};

#endif /* ENZO_ENZO_INITIAL_IMPLOSION2_HPP */

