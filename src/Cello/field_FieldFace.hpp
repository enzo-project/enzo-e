// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldFace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-12
/// @brief    [\ref Field] Interface for the FieldFace class

#ifndef FIELD_FIELD_FACE_HPP
#define FIELD_FIELD_FACE_HPP

class FieldFace {

  /// @class    FieldFace
  /// @ingroup  Field
  /// @brief [\ref Field] Class for loading field faces and storing
  /// field ghosts zones

public: // interface

  /// Constructor
  FieldFace() throw();

#ifdef CONFIG_USE_CHARM
  /// Constructor
  FieldFace(int n, char * array) throw();
#endif

  /// Destructor
  ~FieldFace() throw();

  /// Copy constructor
  FieldFace(const FieldFace & FieldFace) throw();

  /// Assignment operator
  FieldFace & operator= (const FieldFace & FieldFace) throw();

  //----------------------------------------------------------------------

  /// Set whether or not to include ghost zones along each axis
  inline void set_full (bool gx, bool gy = true, bool gz = true)
  {
    full_[axis_x] = gx;
    full_[axis_y] = gy;
    full_[axis_z] = gz;
  }
  
  /// Load from field's face data
  void load(const FieldDescr * field_descr,
	    const FieldBlock * field_block,
	    int                fx,
	    int                fy = 0,
	    int                fz = 0) throw();


  /// Store to field's ghost data
  void store(const FieldDescr * field_descr,
	     FieldBlock *       field_block,
	     int                fx,
	     int                fy = 0,
	     int                fz = 0) throw();

  /// Allocate array_ storage
  void allocate(const FieldDescr * field_descr,
		const FieldBlock * field_block,
		int                fx,
		int                fy = 0,
		int                fz = 0) throw();

  /// Deallocate array_ storage
  void deallocate() throw();

  /// Return the size of the array
  size_t size() const throw() { return array_.size(); };

  /// Return a pointer to the array
  char * array () throw() { return &array_[0]; };

private: // functions

  /// Compute loop limits for load_precision_ and store_precision_
  void loop_limits
  (int *ix0, int *iy0, int *iz0,
   int *nx,  int *ny,  int *nz,
   int nd3[3], int ng3[3],
   int  fx, int  fy, int  fz,
   bool load);


  /// Precision-agnostic function for loading field block face into
  /// the field_face array; returns number of bytes copied
  template<class T>
  size_t load_precision_ (T *       array_face, 
			  const T * field_face,
			  int       nd3[3], 
			  int       ng3[3],
			  int       fx,
			  int       fy = 0,
			  int       fz = 0) throw();

  /// Precision-agnostic function for copying the field_face array into
  /// the field block ghosts; returns number of bytes copied
  template<class T>
  size_t store_precision_ (T *       field_ghosts, 
			   const T * array_ghosts, 
			   int       nd3[3], 
			   int       ng3[3],
			   int       fx,
			   int       fy = 0,
			   int       fz = 0) throw();

private: // attributes

  /// Allocated array used for storing all ghosts and face
  std::vector<char> array_;

  /// Whether to include ghost zones along x,y,z axes
  bool full_[3];

};

#endif /* FIELD_FIELD_FACE_HPP */
