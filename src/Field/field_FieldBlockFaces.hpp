// $Id: field_FieldBlockFaces.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELDBLOCKFACES_HPP
#define FIELD_FIELDBLOCKFACES_HPP

/// @file     field_FieldBlockFaces.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Interface for the FieldBlockFaces class

class FieldBlockFaces {

  /// @class    FieldBlockFaces
  /// @ingroup  Field
  /// @brief    Class for representing and operating on ghost zones

public: // interface

  /// Constructor
  FieldBlockFaces() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~FieldBlockFaces() throw();

  /// Copy constructor
  FieldBlockFaces(const FieldBlockFaces & FieldBlockFaces) throw();

  /// Assignment operator
  FieldBlockFaces & operator= (const FieldBlockFaces & FieldBlockFaces) throw();

  /// Copy ghost zones from FieldBlock to FieldBlockFaces
  void copy_from_block();

  /// Copy ghost zones from FieldBlockFaces to FieldBlock
  void copy_to_block();

  ///  	Initiate a send of ghost zones to another block patch, which
  ///  	may be remote. Nonblocking.
  void send_begin();

  ///  	Complete a send of ghost zones to another block patch, which
  ///  	may be remote. Blocking.
  void send_end();

  ///  	Initiate a receive of ghost zones from another block patch, which
  ///  	may be remote. Nonblocking.
  void recv_begin();

  ///  	Complete a receive of ghost zones from another block patch, which
  ///  	may be remote. Blocking.
  void recv_end();

  ///  	Initiate an exchange of ghost zones associated with a
  ///  	patch. Nonblocking.
  void exchange_begin();

  ///  	Complete an exchange of ghost zones associated with a patch. 
  void exchange_end();

private: // functions


private: // attributes

  /// Pointer to the associated FieldBlock
  FieldBlock * field_block;

  /// Allocated arrays [xm,xp,xm,xp,zm,zp] of ghost values (0 if none)
  char * faces_[6];

  /// Process affinities (e.g. MPI rank, thread id, etc.) of neighbors
  Affinity ** affinity_;

};

#endif /* FIELD_FIELDBLOCKFACES_HPP */

