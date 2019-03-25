#ifndef ENZO_INITIAL_B_CENTER_HPP
#define ENZO_INITIAL_B_CENTER_HPP

// The primary purpose of this class is to initialize the cell-centered B-fields
// after using functions to initialize the face-centered fields the parameter
// file.
//
// Additionally, it provides the static method initialize_bfield_center for
// other Initial classes to call directly.

class EnzoInitialBCenter : public Initial {

  /// @class    EnzoInitialPpmlTest
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializes cell-centered B-field from already
  /// initialized face-centered B-fields

public: // interface

  /// Constructor
  EnzoInitialBCenter (int cycle, double time) throw ();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialBCenter);

  /// CHARM++ migration constructor
  EnzoInitialBCenter(CkMigrateMessage *m)
    : Initial (m)
  {  }


  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    Initial::pup(p);
    // NOTE: update this member function whenever class attributes change
  }

  virtual ~EnzoInitialBCenter() throw()
  {  }

  /// static method that initializes the cell-centered bfield from previously
  /// initialized face-centered bfields
  static void initialize_bfield_center( Block * block );

public: // virtual functions

  /// Initialize a Block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();
};

#endif /* ENZO_INITIAL_B_CENTER_HPP */
