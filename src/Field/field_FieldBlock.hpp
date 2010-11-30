// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_FIELD_BLOCK_HPP
#define FIELD_FIELD_BLOCK_HPP

/// @file     field_FieldBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Oct 12 14:38:21 PDT 2009
/// @todo     Add allocation and deallocation
/// @todo     Implement and test merge(),split()
/// @todo     Implement and test grow(),shrink()
/// @todo     Implement and test read(),write()
/// @todo     Clean allocate_array() and allocate_ghosts() usage
/// @brief    Fortran-style array class.

class FieldBlock {

  /// @class    FieldBlock
  /// @ingroup  Field
  /// @brief Interface between field arrays and low-level (C/fortran)
  /// routines.
  /// 
  /// Defines up to a 4-D fortran-like array for storing 1 or more 3D
  /// arrays.  Axes can be permuted, including the index selecting the
  /// array for storing interleaved arrays.

public: // interface

  /// Create a new uninitialized FieldBlock object
  FieldBlock() throw();

  /// Deconstructor
  ~FieldBlock() throw();

  /// Copy constructor
  FieldBlock(const FieldBlock & field_block) throw ();

  /// Assignment operator
  FieldBlock & operator= (const FieldBlock & field_block) throw ();

  /// Return size of fields on the block, assuming centered
  void size(int * nx, int * ny, int * nz) const throw();

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  void * field_values (int id_field) throw (std::out_of_range);

  /// Return array for the corresponding field, which does not contain
  /// ghosts whether they're allocated or not
  void * field_unknowns (int id_field) throw (std::out_of_range);

  /// Return raw pointer to the array of all fields.  Const since
  /// otherwise dangerous due to varying field sizes, precisions,
  /// padding and alignment
  const char * array ()  const throw () { return array_; };

  /// Return lower and upper+1 index ranges (excluding ghosts)
  void index_range(int * lower_x, int * lower_y, int *lower_z, 
		   int * upper_x, int * upper_y, int *upper_z) const throw ();

  /// Return lower values of the block (excluding ghosts)
  void extent(double * lower_x = 0, double * upper_x = 0, 
		  double * lower_y = 0, double * upper_y = 0,
		  double * lower_z = 0, double * upper_z = 0) const throw ();

  /// Return width of cells along each dimension
  void cell_width(double * hx, double * hy, double * hz) const throw ();

  /// Return the associated field descriptor
  FieldDescr * field_descr() throw ();

  /// Clear specified array(s) to specified value
  void clear(float value = 0.0, 
	     int id_field_first = -1, 
	     int id_field_last  = -1) throw();
 
  /// Return whether array is allocated or not
  bool array_allocated() const throw();

  /// Allocate storage for the field block
  void allocate_array() throw();

  /// Deallocate storage for the field block
  void deallocate_array() throw();

  /// Return whether ghost cells are allocated or not
  bool ghosts_allocated() const throw ();

  /// Allocate and clear ghost values if not already allocated
  void allocate_ghosts() throw ();

  /// Deallocate ghost values if allocated
  void deallocate_ghosts() throw ();

  /// Refresh ghost zones on an internal face
  void refresh_ghosts() throw();

  /// Enforce boundary conditions on a boundary face
  void enforce_boundary(boundary_type boundary, face_type face = face_all) throw();

  /// Split a block into 2, 4, or 8 subblocks; does not delete self
  void split(bool split_x, bool split_y, bool split_z, 
	     FieldBlock ** field_blocks) throw ();

  /// Merge 2, 4, or 8 subblocks into a single block
  FieldBlock * merge(bool merge_x, bool merge_y, bool merge_z, 
		     FieldBlock ** field_blocks) throw ();
 
  /// Read a block from disk.  Create new FieldDescr if not supplied or different.
  /// return NULL iff no new field_descr is created
  FieldDescr * read(File * file, FieldDescr * field_descr = 0) throw ();

  /// Write a block from disk, and optionally associated descriptor
  void write(File * file, FieldDescr * field_descr = 0) const throw ();

  /// Set size of the array block
  void set_size(int nx, int ny=1, int nz=1) throw();

  /// Set array values for a given field
  void set_field_values (int id_field, char * values) throw();

  /// Set the associated field descriptor
  void set_field_descr(FieldDescr * field_descr) throw();

  /// Set the box extent
  void set_extent(double lower_x = 0.0, double upper_x = 1.0, 
		  double lower_y = 0.0, double upper_y = 1.0,
		  double lower_z = 0.0, double upper_z = 1.0) throw();

  //----------------------------------------------------------------------

private: // functions

  /// Given field size, padding, and alignment, compute offset to
  /// start of the next field
  int field_size_adjust_ (int size, int padding, int alignment) const throw();

  /// Given array start and alignment, return first address that is
  /// aligned
  int align_padding_ (char * start, int alignment) const throw()
  { 
    long unsigned start_long = reinterpret_cast<long unsigned>(start);
    return ( alignment - (start_long % alignment) ) % alignment; 
  };

  /// Return the size of the given field, both as (nx,ny,nz) and
  /// return n
  int field_size_ (int id_field, int *nx, int *ny, int *nz) const throw();

  /// Move (not copy) array_ to array and field_values_ to
  /// field_values
  void backup_array_  ( std::vector<char *> & field_values );

  /// Move (not copy) array to array_ and field_values to
  /// field_values_
  void restore_array_ ( std::vector<char *> & field_values )
    throw (std::out_of_range);

  /// Enforce reflecting boundary conditions on a boundary face
  void enforce_boundary_reflecting_(face_type face) throw();
  template<class T>
  void enforce_boundary_reflecting_precision_
  ( face_type face,
    T * array,
    int nx,int ny,int nz,
    int gx,int gy,int gz,
    bool vx,bool vy,bool vz);
  /// Enforce outflow boundary conditions on a boundary face
  void enforce_boundary_outflow_(face_type face) throw();
  /// Enforce inflow boundary conditions on a boundary face
  void enforce_boundary_inflow_(face_type face) throw();
  /// Enforce periodic boundary conditions on a boundary face
  void enforce_boundary_periodic_(face_type face) throw();
  /// Enforce dirichlet boundary conditions on a boundary face
  void enforce_boundary_dirichlet_(face_type face) throw();
  /// Enforce neumann boundary conditions on a boundary face
  void enforce_boundary_neumann_(face_type face) throw();

  //----------------------------------------------------------------------

private: // attributes

  /// Corresponding Field descriptor
  FieldDescr * field_descr_;

  /// Size of fields on the block, assuming centered
  int size_[3];

  /// Allocated array of field values
  char * array_;
  
  /// Pointers into values_ of the first element of each field
  std::vector<char *> field_values_;

  /// Extent of the box associated with the block
  /// WARNING: should not be used for deep AMR due to precision /
  /// range issues
  double lower_[3];
  double upper_[3];

  /// Whether ghost values are allocated or not (make [3] for
  /// directionally split?)
  bool ghosts_allocated_;

};   

#endif /* FIELD_FIELD_BLOCK_HPP */
