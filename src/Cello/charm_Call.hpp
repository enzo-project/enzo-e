// See LICENSE_CELLO file for license and copyright information

/// @file     charm_Call.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-04-08
/// @brief    [\ref Charm] Declaration of the Call class

#ifndef CHARM_CALL_HPP
#define CHARM_CALL_HPP

class Call {

  /// @class    Call
  /// @ingroup  Charm
  /// @brief    [\ref Charm] This class is used for encapsulating a CkCallback
  ///                        along with an associated synchronization.

public: // interface

  /// Constructor
  Call() throw();

  /// Copy constructor
  Call(const Call & call) throw();

  /// Assignment operator
  Call & operator= (const Call & call) throw();

  /// Destructor
  ~Call() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* CHARM_CALL_HPP */

