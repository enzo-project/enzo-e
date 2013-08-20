// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialSedovArray2.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Sat Jun 15 13:43:39 PDT 2013
/// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

#ifndef ENZO_ENZO_INITIAL_SEDOV_ARRAY2_HPP
#define ENZO_ENZO_INITIAL_SEDOV_ARRAY2_HPP

class EnzoInitialSedovArray2 : public Initial {

  /// @class    EnzoInitialSedovArray2
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initialization routine for 2D implosion problem

public: // interface

  /// CHARM++ constructor
  EnzoInitialSedovArray2() throw() { }
  
  /// Constructor
  EnzoInitialSedovArray2(const EnzoConfig * enzo_config) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialSedovArray2);

  /// CHARM++ migration constructor
  EnzoInitialSedovArray2(CkMigrateMessage *m) : Initial (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize the block

  virtual void enforce_block
  ( CommBlock * block, const FieldDescr * field_descr, const Hierarchy * hierarchy ) throw();

private: // attributes

  /// Size of the array of Sedov blasts
  int array_[2];

  /// Relative radius
  double radius_relative_;

  /// Internal and external pressure
  double pressure_in_;
  double pressure_out_;

  /// initial density
  double density_;

  /// Whether PPM or PPML is used
  int hydro_;

};

#endif /* ENZO_ENZO_INITIAL_SEDOV_ARRAY2_HPP */

