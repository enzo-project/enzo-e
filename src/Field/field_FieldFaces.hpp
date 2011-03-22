// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELD_FACES_HPP
#define FIELD_FIELD_FACES_HPP

/// @file     field_FieldFaces.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    [\ref Field] Interface for the FieldFaces class

class FieldFaces {

  /// @class    FieldFaces
  /// @ingroup  Field
  /// @brief    [\ref Field] Class for representing and operating on ghost zones

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

  /// Load from field's face data
  void load(const FieldBlock * field_block,
	    int                field=-1,
	    axis_enum          axis=axis_all, 
	    face_enum          face=face_all) throw();

  /// Copy face data to ghost data
  void copy(const FieldFaces * field_face,
	    int                field=-1,
	    axis_enum          axis=axis_all, 
	    face_enum          face=face_all) throw();

  /// Store to field's ghost data
  void store(FieldBlock *   field_block,
	     int            field=-1,
	     axis_enum      axis=axis_all, 
	     face_enum      face=face_all) throw();

private: // functions

  /// Allocate array_ storage
  void allocate_(FieldBlock * field_block) throw();

  /// Deallocate array_ storage
  void deallocate_() throw();

  /// Index of field, axis, dir into ghosts_ and faces_ arrays
  size_t index_(size_t field, size_t axis, size_t face) 
  { return face + 2*(axis + 3*field); };

  template<class T>
  void load_precision_
  (T * face_values, T * field_values,
   int n[3], int nd[3], int ng[3],
   axis_enum axis, face_enum face );

private: // attributes

  /// Allocated array used for storing all ghosts and faces
  std::vector<char> array_;

  /// Offsets index_[field,axis,face] into array_ of ghost values
  std::vector<int> index_ghosts_;

  /// Offsets index_[field,axis,face] into array_ of faces values
  std::vector<int> index_faces_;

};

#endif /* FIELD_FIELD_FACES_HPP */
