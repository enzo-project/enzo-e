// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldFace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-12
/// @brief    [\ref Field] Interface for the FieldFace class

#ifndef FIELD_FIELD_FACE_HPP
#define FIELD_FIELD_FACE_HPP

class Restrict;
class Prolong;

class FieldFace {

  /// @class    FieldFace
  /// @ingroup  Field
  /// @brief [\ref Field] Class for loading field faces and storing
  /// field ghosts zones

public: // interface

  /// Constructor of uninitialized FieldFace
  FieldFace() throw();

  /// Constructor of initialized FieldFace
  FieldFace (FieldBlock * field_block,
	     const FieldDescr * field_descr) throw();
     
  /// Destructor
  ~FieldFace() throw();

  /// Copy constructor
  FieldFace(const FieldFace & FieldFace) throw();

  /// Assignment operator
  FieldFace & operator= (const FieldFace & FieldFace) throw();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p);

  //----------------------------------------------------------------------

  /// Set whether or not to include ghost zones along each axis
  inline void set_ghost (bool gx, bool gy = true, bool gz = true)
  {
    ghost_[0] = gx;
    ghost_[1] = gy;
    ghost_[2] = gz;
  }

  /// Get whether or not to include ghost zones along each axis
  inline void get_ghost (bool * gx, bool * gy = 0, bool * gz = 0)
  {
    if (gx) (*gx) = ghost_[0];
    if (gy) (*gy) = ghost_[1];
    if (gz) (*gz) = ghost_[2];
  }

  /// Set the face
  inline void set_face (int fx, int fy = 0, int fz = 0)
  {
    face_[0] = fx;
    face_[1] = fy;
    face_[2] = fz;
  }
  
  /// Get the face
  inline void get_face (int * fx, int * fy = 0, int * fz = 0)
  {
    if (fx) (*fx) = face_[0];
    if (fy) (*fy) = face_[1];
    if (fz) (*fz) = face_[2];
  }
  
  /// Create an array with the field's face data
  void load(int * n, char ** array) throw();

  /// Interpolate the data using the given prolongation operator
  void set_prolong(Prolong * prolong, int icx, int icy=0, int icz=0) throw()
  { prolong_ = prolong; 
    set_child_(icx,icy,icz);  }

  /// Restrict the data using the given restriction operator
  void set_restrict(Restrict * restrict, int icx, int icy=0, int icz=0) throw()
  { restrict_ = restrict; 
    set_child_(icx,icy,icz);  }

  /// Copy the input array data to the field's ghost zones
  void store(int n, char * array) throw();

  /// Allocate array_ storage
  char * allocate() throw();

  /// Deallocate array_ storage
  void deallocate() throw();

  /// Return the size of the array
  size_t size() const throw() { return array_.size(); };

  /// Return a pointer to the array
  char * array () throw() 
  { return size() > 0 ?  &array_[0] : 0; };

private: // functions

  /// copy data
  void copy_(const FieldFace & field_face); 

  /// Compute loop limits for store_
  void store_loop_limits_
  (int im3[3], int n3[3], int nd3[3], int ng3[3]);

  /// Compute loop limits for load_
  void load_loop_limits_
  (int im3[3], int n3[3], int nd3[3], int ng3[3]);

  /// Compute loop limits for load_
  void new_loop_limits_
  (int im3[3], int n3[3], int nd3[3], int ng3[3], int op_type);

  void check_new_( int im3[3],int n3[3], int nd3[3], int ng3[3],int op_type);

  /// Set child indices if prolongation or restriction is required
  inline void set_child_ (int icx, int icy = 0, int icz = 0)
  {
    child_[0] = icx;
    child_[1] = icy;
    child_[2] = icz;
  }

  /// Precision-agnostic function for loading field block face into
  /// the field_face array; returns number of bytes copied
  template<class T>
  size_t load_ (T * array_face,  const T * field_face,
			  int nd3[3], int nf3[3], int im3[3]) throw();

  /// Precision-agnostic function for copying the field_face array into
  /// the field block ghosts; returns number of bytes copied
  template<class T>
  size_t store_ (T * field_ghosts,  const T * array_ghosts, 
			   int nd3[3], int nf3[3], int im3[3]) throw();

private: // attributes

  /// field block for this face
  FieldBlock * field_block_;

  /// Field descriptor for the field block
  FieldDescr * field_descr_;

  /// Allocated array used for storing all ghosts and face
  std::vector<char> array_;

  /// Select face, including edges and corners (-1,-1,-1) to (1,1,1)
  int face_[3];

  /// Whether to include ghost zones along x,y,z axes
  bool ghost_[3];

  /// Child index (0,0,0) to (1,1,1) if restriction or prolongation are used
  int child_[3];

  /// Restriction operation if any
  Restrict * restrict_;

  /// Prolongation operation if any
  Prolong * prolong_;

};

#endif /* FIELD_FIELD_FACE_HPP */
