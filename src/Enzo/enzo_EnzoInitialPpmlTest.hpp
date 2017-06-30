// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialPpmlTest.hpp
/// @author   First Last (email@address)
/// @date     2017-06-12
/// @brief    [\ref Enzo] Declaration of the EnzoInitialPpmlTest class

#ifndef ENZO_PPML_TEST_HPP
#define ENZO_PPML_TEST_HPP

class EnzoInitialPpmlTest : public Initial {

  /// @class    EnzoInitialPpmlTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the PpmlTest problem

public: // interface

  /// Constructor
  EnzoInitialPpmlTest (int cycle, double time,
		       const EnzoConfig * enzo_config) throw ();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialPpmlTest);

    /// CHARM++ migration constructor
  EnzoInitialPpmlTest(CkMigrateMessage *m)
    : Initial (m)
  {  }


  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    Initial::pup(p);

    // NOTE: update this member function whenever class attributes change
  }

  ~EnzoInitialPpmlTest() throw()
  {  } 

public: // virtual functions

  /// Initialize a Block
  virtual void enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw();

private: // functions

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* ENZO_PPML_TEST_HPP */

