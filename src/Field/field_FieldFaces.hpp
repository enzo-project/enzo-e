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

private:

  /// FieldFaces must be allocated with field_block argument

  FieldFaces() throw();

public: // interface

  /// Constructor
  FieldFaces(FieldBlock * field_block) throw();

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

  /// Copy in face zones from a field block
  void load() throw();

  /// Copy out ghost zones to a field block
  void save() throw();


private: // functions

  /// Allocate array_ storage

  void allocate_() throw();

  /// Deallocate array_ storage

  void deallocate_() throw();

private: // attributes

  /// Pointer to the corresponding FieldBlock
  FieldBlock * field_block_;

  /// Arrays [field][axis][dir] of ghost values
  std::vector<char *> ghosts_;

  /// Arrays [field][axis][dir] of inner-face values
  std::vector<char *> faces_;

  /// Allocated array used for storing all ghosts and faces
  char * array_;

};

#endif /* FIELD_FIELD_FACES_HPP */
