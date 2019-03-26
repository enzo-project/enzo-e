// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialFeedbackTest.hpp
/// @author   Andrew Emerick (emerick@astro.columbia.edu)
/// @date
/// @brief    [\ref Enzo] Initialization routine for Feedback test problem

#ifndef ENZO_ENZO_INITIAL_FEEDBACK_TEST_HPP
#define ENZO_ENZO_INITIAL_FEEDBACK_TEST_HPP

class EnzoInitialFeedbackTest : public Initial {

  /// @class     EnzoInitialFeedbackTest
  /// @ingroup   Enzo
  /// @brif      [\ref Enzo] Initialization routine for feedback test problem

public:  // interface

  ///Charm++ constructor
  EnzoInitialFeedbackTest(const EnzoConfig * enzo_config) throw();

  /// Charm++ PUP::able declaration
  PUPable_decl(EnzoInitialFeedbackTest);

  /// CHARM++ migration constructor
  EnzoInitialFeedbackTest(CkMigrateMessage *m)
      : Initial (m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  // Initialize the block

  virtual void enforce_block
  (   Block * block, const Hierarchy * hierarchy) throw();

  // Destructor
  virtual ~EnzoInitialFeedbackTest() throw() {};

private:

};

#endif /* ENZO_ENZO_INITIAL_FEEDBACK_TEST_HPP */
