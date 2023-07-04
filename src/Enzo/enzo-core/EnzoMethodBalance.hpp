// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodBalance.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     <2022-12-12 Mon>
/// @brief    [\ref Enzo] Declaration for the EnzoMethodBalance class

#ifndef ENZO_ENZO_METHOD_BALANCE_HPP
#define ENZO_ENZO_METHOD_BALANCE_HPP

class EnzoMethodBalance : public Method {

  /// @class    EnzoMethodBalance
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's PPM hydro method

public: // interface

  /// Create a new EnzoMethodBalance object
  EnzoMethodBalance();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodBalance);

  /// Charm++ PUP::able migration constructor
  EnzoMethodBalance (CkMigrateMessage *m)
    : Method (m), ip_next_(-1)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void do_migrate(EnzoBlock * enzo_block);
  void done(EnzoBlock * enzo_block);

public: // virtual methods

  /// Apply the method to advance a block one timestep 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "balance"; }

protected: // attributes

  /// Process to migrate to
  int ip_next_;

};

#endif /* ENZO_ENZO_METHOD_BALANCE_HPP */
