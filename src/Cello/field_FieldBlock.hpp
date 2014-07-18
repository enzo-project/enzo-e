// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Oct 12 14:38:21 PDT 2009
/// @brief    [\ref Field] Fortran-style array class.

#ifndef FIELD_FIELD_BLOCK_HPP
#define FIELD_FIELD_BLOCK_HPP

class CommBlock;
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
  FieldBlock(FieldDescr * field_descr = 0,
	     int nx=0, int ny=1, int nz=1) throw();

  /// Deconstructor
  ~FieldBlock() throw();

  /// Copy constructor
  FieldBlock(const FieldBlock & field_block) throw ();

  /// Assignment operator
  FieldBlock & operator= (const FieldBlock & field_block) throw ();

  void pup(PUP::er &p) ;

  /// Return size of fields on the block, assuming centered
  void size(int * nx, int * ny = 0, int * nz = 0) const throw();

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  char * field_values (int id_field) throw (std::out_of_range);
  char * field_values (std::string name) throw (std::out_of_range)
  { return field_values (field_descr_->field_id(name)); }

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  const char * field_values (int id_field) const throw (std::out_of_range);
  const char * field_values (std::string name) const throw (std::out_of_range)
  { return field_values (field_descr_->field_id(name)); }

  /// Return array for the corresponding field, which does not contain
  /// ghosts whether they're allocated or not
  char * field_unknowns ( int id_field) throw (std::out_of_range);
  char * field_unknowns (std::string name) throw (std::out_of_range)
  { return field_unknowns (field_descr_->field_id(name)); }

  const char * field_unknowns ( int id_field) const throw (std::out_of_range);
  const char * field_unknowns (std::string name) const throw (std::out_of_range)
  { return field_unknowns (field_descr_->field_id(name)); }

  /// Return raw pointer to the array of all fields.  Const since
  /// otherwise dangerous due to varying field sizes, precisions,
  /// padding and alignment
  const char * array ()  const throw () 
  { return array_allocated() ? &array_[0] : NULL; };

  /// Return width of cells along each dimension
  void cell_width(double xm,   double xp,   double * hx,
		  double ym=0, double yp=0, double * hy=0,
		  double zm=0, double zp=0, double * hz=0) const throw ();

  /// Clear specified array(s) to specified value
  void clear ( float value = 0.0, 
	       int id_field_first = -1, 
	       int id_field_last  = -1) throw();
 
  /// Return whether array is allocated or not
  bool array_allocated() const throw()
  { return array_.size() > 0; }

  /// Return whether array is allocated or not
  size_t array_size() const throw()
  { return array_.size(); }

  /// Allocate storage for the field block
  void allocate_array(bool ghosts_allocated = false) throw();

  /// Reallocate storage for the field block, e.g. when changing
  /// from ghosts to non-ghosts [ costly for large blocks ]
  void reallocate_array(bool ghosts_allocated = false) throw();

  /// Deallocate storage for the field block
  void deallocate_array() throw();

  /// Return whether ghost cells are allocated or not.  
  bool ghosts_allocated() const throw ()
  {  return ghosts_allocated_; }

  /// Return the number of elements (nx,ny,nz) along each axis, and total
  /// number of bytes n
  int field_size (int id_field, int *nx=0, int *ny=0, int *nz=0) const throw();

  /// Multiply FieldBlock 1 by FieldBlock 2
  void mul (int id_1, int id_2);
  /// Divide FieldBlock 1 by FieldBlock 2
  void div (int id_1, int id_2);
  /// Multiply FieldBlock by a constant
  void scale (int id_1, double scale);


  const FieldDescr * field_descr() { return field_descr_; }

  //----------------------------------------------------------------------

  /// Print basic field characteristics for debugging
  void print (const char * message,
	      bool use_file = false) const throw();

  //----------------------------------------------------------------------
  // Functions passed on to FieldDescr

  /// Return the number of fields
  int field_count() const throw() 
  { return field_descr_->field_count(); }

  /// Return name of the ith field
  std::string field_name(size_t id_field) const throw(std::out_of_range)
  { return field_descr_->field_name(id_field); }

  /// Return whether the field has been inserted
  bool is_field(const std::string & name) const throw()
  { return field_descr_->is_field (name); }

  /// Return the integer handle for the named field
  int field_id(const std::string & name) const throw()
  { return field_descr_->field_id(name); }

  //----------------------------------------------------------------------
  // Groups
  //----------------------------------------------------------------------

  //----------------------------------------------------------------------
  // Groups
  //----------------------------------------------------------------------

  /// Add the field to a group
  void add_to_group(const std::string field,
		    const std::string group) throw(std::out_of_range);

  /// Return whether the given field is in the given group
  bool is_in_group(const std::string field,
		   const std::string group) const throw(std::out_of_range);

  //----------------------------------------------------------------------
  // Properties
  //----------------------------------------------------------------------

  /// alignment in bytes of fields in memory
  int alignment() const throw()
  { return field_descr_->alignment(); }

  /// padding in bytes between fields in memory
  int padding() const throw()
  { return field_descr_->padding(); }

  /// courant number for fields
  double courant() const throw()
  { return field_descr_->courant(); }

  /// centering of given field
  void centering(int id_field, bool * cx, bool * cy = 0, bool * cz = 0) const 
    throw(std::out_of_range)
  { return field_descr_->centering(id_field,cx,cy,cz); }

  /// depth of ghost zones of given field
  void ghosts(int id_field, int * gx, int * gy = 0, int * gz = 0) const 
    throw(std::out_of_range)
  { return field_descr_->ghosts(id_field,gx,gy,gz); }

  /// Return precision of given field
  precision_type precision(int id_field) const throw(std::out_of_range)
  { return field_descr_->precision(id_field); }

  /// Number of bytes per element required by the given field
  int bytes_per_element(int id_field) const throw()
  { return field_descr_->bytes_per_element(id_field); }

  /// minimum value for the field
  double minimum_value(int id_field) const throw(std::out_of_range)
  { return field_descr_->minimum_value(id_field); }

  /// minimum action for the field
  field_action_type minimum_action(int id_field) const
    throw(std::out_of_range)
  { return field_descr_->minimum_action(id_field); }

  /// maximum value for the field
  double maximum_value(int id_field) const  throw(std::out_of_range)
  { return field_descr_->maximum_value(id_field); }

  /// maximum action for the field
  field_action_type maximum_action(int id_field) const 
    throw(std::out_of_range)
  { return field_descr_->maximum_action(id_field); }

  //--------------------------------------------------
private: // functions
  //--------------------------------------------------

  /// Given field size and padding, compute offset to start of the next field
  int adjust_padding_ (int size, int padding) const throw();

  /// Given field size and alignment, compute offset to start of the next field
  int adjust_alignment_ (int size, int alignment) const throw();

  /// Given array start and alignment, return first address that is
  /// aligned
  int align_padding_ (int alignment) const throw();

  /// Move (not copy) array to array_ and offsets to
  /// offsets_
  void restore_array_ ( const char       * array_from,
			std::vector<int> & offsets_from )
    throw (std::out_of_range);

  template <class T>
  void print_
  (const T * field,
   const char * field_name,
   const char * message,
   //   double lower [3],
   FILE * fp,
   int ixm,int iym,int izm,
   int ixp,int iyp,int izp,
   int nx, int ny, int nz,
   int gx, int gy ,int gz,
   //   double hx, double hy ,double hz,
   int nxd,int nyd) const;

  template<class T>
  void mul_ (T * field_1, T * field_2);

  template<class T>
  void div_ (T * field_1, T * field_2);

  template<class T>
  void scale_ (T * field, double value);

  //----------------------------------------------------------------------

  /// Set alignment
  void set_alignment(int alignment) throw();

  /// Set padding
  void set_padding(int padding) throw();

  /// Set courant
  void set_courant(double courant) throw();

  /// Set centering for a field
  void set_centering(int id_field, bool cx, bool cy=true, bool cz=true) 
    throw(std::out_of_range);

  /// Set ghosts for a field
  void set_ghosts(int id_field, int gx, int gy=0, int gz=0) 
    throw(std::out_of_range);

  /// Set precision for a field
  void set_precision(int id_field, precision_type precision) 
    throw(std::out_of_range);

  /// Set minimum bound and action
  void set_minimum (int id_field, double min_value,
		    field_action_type min_action) 
    throw(std::out_of_range);

  /// Set maximum bound and action
  void set_maximum (int id_field, double max_value, 
		    field_action_type max_action) 
    throw(std::out_of_range);

  /// Insert a new field
  int insert_field(const std::string & name_field) throw();

private: // attributes

  /// Pointer to the associated field descriptor
  FieldDescr * field_descr_;

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
