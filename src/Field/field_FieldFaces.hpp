// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELD_FACES_HPP
#define FIELD_FIELD_FACES_HPP

/// @file     field_FieldFaces.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Interface for the FieldFaces class

class FieldFaces {

  /// @class    FieldFaces
  /// @ingroup  Field
  /// @brief    Class for representing and operating on ghost zones

public: // interface

  /// Constructor
  FieldFaces() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~FieldFaces() throw();

  /// Copy constructor
  FieldFaces(const FieldFaces & FieldFaces) throw();

  /// Assignment operator
  FieldFaces & operator= (const FieldFaces & FieldFaces) throw();

  //----------------------------------------------------------------------

  /// Copy ghost zones from FieldBlock to FieldFaces
  void copy_from_block();

  /// Copy ghost zones from FieldFaces to FieldBlock
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

  /// Size of the arrays
  int n_;

  /// Allocated arrays [xm,xp,xm,xp,zm,zp] of ghost values (0 if none)
  char * faces_[6];

};

#endif /* FIELD_FIELD_FACES_HPP */
