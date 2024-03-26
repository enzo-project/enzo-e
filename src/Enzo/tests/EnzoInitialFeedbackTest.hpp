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
  EnzoInitialFeedbackTest(int cycle, double time, ParameterGroup p) throw();

  /// Charm++ PUP::able declaration
  PUPable_decl(EnzoInitialFeedbackTest);

  /// CHARM++ migration constructor
  EnzoInitialFeedbackTest(CkMigrateMessage *m)
    : Initial (m),
      density_CGS_(0.0),
      temperature_(0.0),
      metal_fraction_(0.0),
      species_densities_CGS_(),
      num_particles(0),
      position{},
      mass(),
      luminosity()
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  // Initialize the block

  virtual void enforce_block
  (   Block * block, const Hierarchy * hierarchy) throw();

  // Destructor
  virtual ~EnzoInitialFeedbackTest() throw() {};

private:
  
  double density_CGS_;
  double temperature_;
  double metal_fraction_;
  std::map<std::string, double> species_densities_CGS_;
  int num_particles;
  std::vector<double> position[3];
  std::vector<double> mass;
  std::vector<double> luminosity;

};

#endif /* ENZO_ENZO_INITIAL_FEEDBACK_TEST_HPP */
