// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialZeldovichPancake.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 22 2022
/// @brief    [\ref Enzo] Initialization routine ZeldovichPancake cosmology
///           test problem.

#ifndef ENZO_ENZO_INITIAL_ZELDOVICH_PANCAKE_HPP
#define ENZO_ENZO_INITIAL_ZELDOVICH_PANCAKE_HPP

class EnzoInitialZeldovichPancake : public Initial {
  /// @class    EnzoInitialZeldovichPancake
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for the axis-aligned Zeldovich Pancake
  /// test problem.
  ///
  /// In the future, we might want to combine with the EnzoInitialInclinedWave
  /// initializer or reuse some of the machinery so that we can incline this
  /// problem.

public: // interface

  /// Constructor
  EnzoInitialZeldovichPancake(int cycle, double time,
                              std::string aligned_ax_name);

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialZeldovichPancake);

  /// CHARM++ migration constructor
  EnzoInitialZeldovichPancake(CkMigrateMessage *m)
    : Initial (m), aligned_ax_(0)
  {  }

  /// Destructor
  virtual ~EnzoInitialZeldovichPancake() throw()
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) {
    // NOTE: update whenever attributes change
    TRACEPUP;
    Initial::pup(p);
    p | aligned_ax_;
  }

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // attributes

  /// specify orientation of problem 0,1,2 correspond to x,y,z
  int aligned_ax_;
};

#endif /* ENZO_ENZO_INITIAL_ZELDOVICH_PANCAKE_HPP */
