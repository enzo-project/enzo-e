// $Id: field_FieldGhosts.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELDGHOSTS_HPP
#define FIELD_FIELDGHOSTS_HPP

/// @file     field_FieldGhosts.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Interface for the FieldGhosts class

class FieldGhosts {

  /// @class    FieldGhosts
  /// @ingroup  Field
  /// @brief    Class for representing and operating on ghost zones

public: // interface

  /// Constructor
  FieldGhosts() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~FieldGhosts() throw();

  /// Copy constructor
  FieldGhosts(const FieldGhosts & FieldGhosts) throw();

  /// Assignment operator
  FieldGhosts & operator= (const FieldGhosts & FieldGhosts) throw();

  /// Copy ghost zones from FieldBlock to FieldGhosts
  void copy_from_block();

  /// Copy ghost zones from FieldGhosts to FieldBlock
  void copy_to_block();

  ///  	Initiate a send of ghost zones to another block patch, which
  ///  	may be remote. Nonblocking.
  void send();

  ///  	Complete a send of ghost zones to another block patch, which
  ///  	may be remote. Blocking.
  void send_wait();

  ///  	Initiate a receive of ghost zones from another block patch, which
  ///  	may be remote. Nonblocking.
  void recv();

  ///  	Complete a receive of ghost zones from another block patch, which
  ///  	may be remote. Blocking.
  void recv_wait();

  ///  	Initiate an exchange of ghost zones associated with a
  ///  	patch. Nonblocking.
  void exchange();

  ///  	Complete an exchange of ghost zones associated with a patch. 
  void exchange_wait();

private: // functions


private: // attributes

  /// Pointer to the associated FieldBlock
  FieldBlock * field_block;

  /// Allocated arrays [axis][face] of ghost values
  char *** ghosts_;

  /// Process affinities (e.g. MPI rank, thread id, etc.) of neighbors
  Affinity ** affinity_;

};

#endif /* FIELD_FIELDGHOSTS_HPP */

