// $Id: field_FieldFace.hpp 2179 2011-04-06 22:40:24Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELD_FACE_HPP
#define FIELD_FIELD_FACE_HPP

/// @file     field_FieldFace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    [\ref Field] Interface for the FieldFace class

class FieldFace {

  /// @class    FieldFace
  /// @ingroup  Field
  /// @brief [\ref Field] Class for loading field faces and storing
  /// field ghosts zones

public: // interface

  /// Constructor
  FieldFace() throw();

  /// Constructor
  FieldFace(int n, char * array) throw()
    : array_()
  {
    TRACE("FieldFace");
    array_.resize(n);
    TRACE("FieldFace");
    printf("n=%d %g\n",n,*((float *)(&array_[0])));
    for (int i=0; i<n; i++) array_[i] = array[i]; 
    TRACE("FieldFace");
  };

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~FieldFace() throw();

  /// Copy constructor
  FieldFace(const FieldFace & FieldFace) throw();

  /// Assignment operator
  FieldFace & operator= (const FieldFace & FieldFace) throw();

  //----------------------------------------------------------------------

  /// Load from field's face data
  void load(const FieldDescr * field_descr,
	    const FieldBlock * field_block,
	    axis_enum          axis,
	    face_enum          face,
	    bool               full) throw();


  /// Store to field's ghost data
  void store(const FieldDescr * field_descr,
	     FieldBlock *       field_block,
	     axis_enum          axis,
	     face_enum          face,
	     bool               full) throw();

  /// Return the size of the array
  size_t size() const throw() { return array_.size(); };

  /// Return a pointer to the array
  char * array () throw() { return &array_[0]; };

  /// Print basic field characteristics for debugging
  void print (const FieldDescr * field_descr,
	      const FieldBlock * field_block,
	      axis_enum axis, face_enum face,
	      const char * message = 0) const throw();

private: // functions

  /// Allocate array_ storage
  void allocate_(const FieldDescr * field_descr,
		 const FieldBlock * field_block,
		 axis_enum          axis,
		 bool               full) throw();

  /// Deallocate array_ storage
  void deallocate_() throw();

  /// Precision-agnostic function for loading field block face into
  /// the field_face array; returns number of bytes copied
  template<class T>
  size_t load_precision_
  (T * array_face, const T * field_face,
   int n[3], int nd[3], int ng[3],
   axis_enum axis, face_enum face, bool full ) throw();

  /// Precision-agnostic function for copying the field_face array into
  /// the field block ghosts; returns number of bytes copied
  template<class T>
  size_t store_precision_
  (T * field_ghosts, const T * array_ghosts, 
   int n[3], int nd[3], int ng[3],
   axis_enum axis, face_enum face, bool full ) throw();

private: // attributes

  /// Allocated array used for storing all ghosts and face
  std::vector<char> array_;

};

#endif /* FIELD_FIELD_FACE_HPP */
