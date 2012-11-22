// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Oct 12 14:38:21 PDT 2009
/// @brief    [\ref Field] Fortran-style array class.

#ifndef FIELD_FIELD_BLOCK_HPP
#define FIELD_FIELD_BLOCK_HPP

class Patch;
class Block;
class FieldDescr;

class FieldBlock {

  /// @class    FieldBlock
  /// @ingroup  Field
  /// @brief [\ref Field] Interface between field arrays and low-level
  /// (C/fortran) routines.
  /// 
  /// Defines up to a 4-D fortran-like array for storing 1 or more 3D
  /// arrays.  Axes can be permuted, including the index selecting the
  /// array for storing interleaved arrays.

  friend class FieldFace;

public: // interface

  /// Create a new uninitialized FieldBlock object
  FieldBlock(int nx=0, int ny=1, int nz=1) throw();

  /// Deconstructor
  ~FieldBlock() throw();

  /// Copy constructor
  FieldBlock(const FieldBlock & field_block) throw ();

  /// Assignment operator
  FieldBlock & operator= (const FieldBlock & field_block) throw ();

#ifdef CONFIG_USE_CHARM

  void pup(PUP::er &p) 
  {

    TRACEPUP;

    TRACE1("this = %p",this);
    PUParray(p,size_,3);

    p | array_;
    p | offsets_;
    p | ghosts_allocated_;

    TRACE3("size = %d %d %d\n",size_[0],size_[1],size_[2]);
    TRACE1("array_.size() = %d",array_.size());
    TRACE1("offsets_.size() = %d",offsets_.size());
  }

#endif

  /// Return size of fields on the block, assuming centered
  void size(int * nx = 0, int * ny = 0, int * nz = 0) const throw();

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  char * field_values (int id_field) throw (std::out_of_range);

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  const char * field_values (int id_field) const throw (std::out_of_range);

  /// Return array for the corresponding field, which does not contain
  /// ghosts whether they're allocated or not
  char * field_unknowns ( const FieldDescr * field_descr,
			  int id_field) throw (std::out_of_range);

  const char * field_unknowns ( const FieldDescr * field_descr,
				int id_field) const throw (std::out_of_range);

  /// Return raw pointer to the array of all fields.  Const since
  /// otherwise dangerous due to varying field sizes, precisions,
  /// padding and alignment
  const char * array ()  const throw () 
  { return array_allocated() ? &array_[0] : NULL; };

  /// Return width of cells along each dimension
  void cell_width(double xm,   double xp,   double * hx,
		  double ym=0, double yp=0, double * hy=0,
		  double zm=0, double zp=0, double * hz=0) const throw ();

  // /// Return the associated field faces object
  // FieldFaces * field_faces(const FieldDescr * field_descr) throw ();

  /// Clear specified array(s) to specified value
  void clear ( const FieldDescr * field_descr,
	       float value = 0.0, 
	       int id_field_first = -1, 
	       int id_field_last  = -1) throw();
 
  /// Return whether array is allocated or not
  bool array_allocated() const throw();

  /// Allocate storage for the field block
  void allocate_array(const FieldDescr * field_descr) throw();

  /// Deallocate storage for the field block
  void deallocate_array() throw();

  /// Return whether ghost cells are allocated or not
  bool ghosts_allocated() const throw ();

  /// Set whether allocate will allocate array with ghosts or not
  void set_ghosts_allocated(bool lg) throw ()
  { ghosts_allocated_ = lg; }

  // /// Allocate and clear ghost values if not already allocated
  // void allocate_ghosts(const FieldDescr * field_descr) throw ();

  // /// Deallocate ghost values if allocated
  // void deallocate_ghosts(const FieldDescr * field_descr) throw ();


#ifndef CONFIG_USE_CHARM /* MPI only */

  /// Refresh ghost zones on an internal face
  void refresh_ghosts(const FieldDescr * field_descr,
		      const GroupProcess * group_process,
		      const Layout * layout,
		      int ibx, int iby, int ibz,
		      int fx,  int fy,  int fz) throw();

#endif

  /// Split a block into 2, 4, or 8 subblocks; does not delete self
  void split(bool split_x, bool split_y, bool split_z, 
	     FieldBlock ** field_blocks) throw ();

  /// Merge 2, 4, or 8 subblocks into a single block
  FieldBlock * merge(bool merge_x, bool merge_y, bool merge_z, 
		     FieldBlock ** field_blocks) throw ();
 
  /// Return the number of elements (nx,ny,nz) along each axis, and total
  /// number of bytes n
  int field_size (const FieldDescr * field_descr, int id_field, 
		  int *nx = 0, int *ny = 0, int *nz = 0) const throw();

  //----------------------------------------------------------------------

  /// Print basic field characteristics for debugging
  void print (const FieldDescr * field_descr,
	      const char * message,
	      double lower[3], double upper[3],
	      bool use_file = false) const throw();

private: // functions

  /// Given field size and padding, compute offset to start of the next field
  int adjust_padding_ (int size, int padding) const throw();

  /// Given field size and alignment, compute offset to start of the next field
  int adjust_alignment_ (int size, int alignment) const throw();

  /// Given array start and alignment, return first address that is
  /// aligned
  int align_padding_ (int alignment) const throw();

  /// Move (not copy) array_ to array and offsets_ to
  /// offsets
  void backup_array_  ( const FieldDescr * field_descr,
			std::vector<int> & offsets );

  /// Move (not copy) array to array_ and offsets to
  /// offsets_
  void restore_array_ ( const FieldDescr * field_descr,
			const char       * array_from,
			std::vector<int> & offsets_from )
    throw (std::out_of_range);

  //----------------------------------------------------------------------

private: // attributes

  /// Size of fields on the block, assuming centered
  int size_[3];

  /// Allocated array of field values
  std::vector<char> array_;

  /// Offsets into values_ of the first element of each field
  std::vector<int> offsets_;

  /// Whether ghost values are allocated or not 
  bool ghosts_allocated_;

};   

#endif /* FIELD_FIELD_BLOCK_HPP */
