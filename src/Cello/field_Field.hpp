// See LICENSE_CELLO file for license and copyright information

/// @file     field_Field.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Sep  8 14:16:32 PDT 2014
/// @brief    [\ref Field] Declaration of the Field class
///
/// The Field class is used to unify the interface of the global FieldDescr
/// object and a given FieldBlock object.

#ifndef FIELD_FIELD_HPP
#define FIELD_FIELD_HPP

class Field {

  /// @class    Field
  /// @ingroup  Field
  /// @brief    [\ref Field] 

public: // interface

  /// Constructor
  Field(FieldDescr * field_descr,
	FieldBlock * field_block) throw()
    : field_descr_ (field_descr),
      field_block_ (field_block)
  {}

  /// Copy constructor
  Field(const Field & field) throw()
  {
    field_descr_ = field.field_descr_;
    field_block_ = field.field_block_; 
  }

  /// Assignment operator
  Field & operator= (const Field & field) throw()
  { 
    field_descr_ = field.field_descr_;
    field_block_ = field.field_block_;
    return *this;
  }

  /// Destructor
  ~Field() throw()
  {};

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    WARNING ("Field::pup()",
	     "Skipping since Field is intended as transient objects");
  };
  
  /// Return the field descriptor for this field
  FieldDescr * field_descr() { return field_descr_; }

  /// Return the field block for this field
  FieldBlock * field_block() { return field_block_; }

  //==================================================
  // FieldDescr
  //==================================================

  /// Set alignment
  void set_alignment(int alignment) throw()
  { field_descr_->set_alignment(alignment); }

  /// Set padding
  void set_padding(int padding) throw()
  { field_descr_->set_padding(padding); }

  /// Set courant
  void set_courant(double courant) throw()
  { field_descr_->set_courant(courant); }

  /// Set centering for a field
  void set_centering(int id, int cx, int cy=0, int cz=0) 
    throw(std::out_of_range)
  { field_descr_->set_centering(id,cx,cy,cz); }

  /// Set ghosts for a field
  void set_ghosts(int id, int gx, int gy=0, int gz=0) 
    throw(std::out_of_range)
  { field_descr_->set_ghosts(id,gx,gy,gz); }


  /// Set precision for a field
  void set_precision(int id, precision_type precision) 
    throw(std::out_of_range)
  { field_descr_->set_precision(id,precision); }

  /// Insert a new field
  int insert_field(const std::string & name) throw()
  { return field_descr_->insert_field(name); }

  /// Return the number of fields
  int field_count() const throw()
  { return field_descr_->field_count(); }

  /// Return name of the ith field
  std::string field_name(size_t id) const throw(std::out_of_range)
  { return field_descr_->field_name(id); }

  /// Return whether the field has been inserted
  bool is_field(const std::string & name) const throw()
  { return field_descr_->is_field(name); }

  /// Return the integer handle for the named field
  int field_id(const std::string & name) const throw()
  { return field_descr_->field_id(name); }

  //----------------------------------------------------------------------
  // Properties
  //----------------------------------------------------------------------

  Grouping * groups () 
  { return field_descr_->groups(); }

  const Grouping * groups () const 
  { return field_descr_->groups(); }

  /// alignment in bytes of fields in memory
  int alignment() const throw()
  { return field_descr_-> alignment();}

  /// padding in bytes between fields in memory
  int padding() const throw()
  { return field_descr_->padding() ;}

  /// courant number for fields
  double courant() const throw()
  { return field_descr_->courant() ;}

  /// centering of given field
  void centering(int id, int * cx, int * cy = 0, int * cz = 0) const 
    throw(std::out_of_range)
  { return field_descr_->centering(id,cx,cy,cz); }

  /// depth of ghost zones of given field
  void ghosts(int id, int * gx, int * gy = 0, int * gz = 0) const 
    throw(std::out_of_range)
  { return field_descr_->ghosts(id,gx,gy,gz); }

  /// Return precision of given field
  int precision(int id) const throw(std::out_of_range)
  { return field_descr_->precision(id); }

  /// Number of bytes per element required by the given field
  int bytes_per_element(int id) const throw()
  { return field_descr_->bytes_per_element(id); }

  //==================================================
  // FieldBlock
  //==================================================

  /// Return size of fields on the block, assuming centered
  void size(int * nx, int * ny = 0, int * nz = 0) const throw()
  { field_block_->size(nx,ny,nz); }

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  char * values (int id) throw (std::out_of_range)
  { return field_block_->values(id); }

  char * values (std::string name) throw (std::out_of_range)
  { return field_block_->values(name); }

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  const char * values (int id) const throw (std::out_of_range)
  { return field_block_->values(id); }

  const char * values (std::string name) const throw (std::out_of_range)
  { return field_block_->values(name); }

  /// Return array for the corresponding field, which does not contain
  /// ghosts whether they're allocated or not
  char * unknowns ( int id) throw (std::out_of_range)
  { return field_block_->unknowns(id); }

  char * unknowns (std::string name) throw (std::out_of_range)
  { return field_block_->unknowns(name); }

  const char * unknowns ( int id) const throw (std::out_of_range)
  { return field_block_->unknowns(id); }

  const char * unknowns (std::string name) const throw (std::out_of_range)
  { return field_block_->unknowns(name); }

  /// Return raw pointer to the array of all fields.  Const since
  /// otherwise dangerous due to varying field sizes, precisions,
  /// padding and alignment
  const char * array ()  const throw () 
  { return field_block_->array(); }

  /// Return width of cells along each dimension
  void cell_width(double xm,   double xp,   double * hx,
		  double ym=0, double yp=0, double * hy=0,
		  double zm=0, double zp=0, double * hz=0) const throw ()
  { field_block_->cell_width(xm,xp,hx,
			     ym,yp,hy,
			     zm,zp,hz); }

  /// Clear specified array(s) to specified value
  void clear ( float value = 0.0, 
	       int id_first = -1, 
	       int id_last  = -1) throw()
  { field_block_->clear (value,id_first,id_last); }
 
  /// Return whether array is allocated or not
  bool array_allocated() const throw()
  { return field_block_->array_allocated(); }

  /// Return whether array is allocated or not
  size_t array_size() const throw()
  { return field_block_->array_size(); }

  /// Allocate storage for the field block
  void allocate_array(bool ghosts_allocated = false) throw()
  { field_block_->allocate_array(ghosts_allocated); }

  /// Reallocate storage for the field block, e.g. when changing
  /// from ghosts to non-ghosts [ costly for large blocks ]
  void reallocate_array(bool ghosts_allocated = false) throw()
  { field_block_->reallocate_array(ghosts_allocated); }

  /// Deallocate storage for the field block
  void deallocate_array() throw()
  { field_block_->deallocate_array(); }

  /// Return whether ghost cells are allocated or not.  
  bool ghosts_allocated() const throw ()
  { return field_block_->ghosts_allocated(); }

  /// Return the number of elements (nx,ny,nz) along each axis, and total
  /// number of bytes n
  int field_size (int id, int *nx=0, int *ny=0, int *nz=0) const throw()
  { return field_block_->field_size(id,nx,ny,nz); }

  //----------------------------------------------------------------------

  /// Print basic field characteristics for debugging
  void print (const char * message,
	      bool use_file = false) const throw()
  { field_block_->print(message,use_file); }

private: // attributes

  /// Field descriptor for global field data
  FieldDescr * field_descr_;

  /// Field block for specific field data
  FieldBlock * field_block_;

  // NOTE: change pup() function whenever attributes change

};

#endif /* FIELD_FIELD_HPP */

