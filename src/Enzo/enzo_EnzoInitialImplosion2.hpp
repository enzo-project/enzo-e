// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialImplosion2.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:26:38 PST 2011
/// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

#ifndef ENZO_ENZO_INITIAL_IMPLOSION2_HPP
#define ENZO_ENZO_INITIAL_IMPLOSION2_HPP

class EnzoInitialImplosion2 : public Initial {

  /// @class    EnzoInitialImplosion2
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

public: // interface

  /// Constructor
  EnzoInitialImplosion2(int cycle, double time) throw();

  /// Initialize the block

  virtual void enforce (Hierarchy * hierarchy,
			const FieldDescr * field_descr,
			Block * block = NULL) throw();

};

#endif /* ENZO_ENZO_INITIAL_IMPLOSION2_HPP */

