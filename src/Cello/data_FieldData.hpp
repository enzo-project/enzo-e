// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Oct 12 14:38:21 PDT 2009
/// @brief    [\ref Data] Fortran-style array class.

#ifndef DATA_FIELD_DATA_HPP
#define DATA_FIELD_DATA_HPP

class Block;
class FieldDescr;

class FieldData {

  /// @class    FieldData
  /// @ingroup  Data
  /// @brief [\ref Data] Interface between field arrays and low-level
  /// (C/fortran) routines.
  /// 
  /// A FieldData stores up to a 4D fortran-like array for
  /// permanently storing 1 or more 3D arrays.  Axes can be permuted,
  /// including the index selecting the array for storing interleaved
  /// arrays.  Temporary fields can be allocated and deallocated as
  /// well, and are stored in a vector of pointers to 3D arrays.  Descriptive
  /// information, such as number of ghost zones, padding, centering, etc.,
  /// are stored in a separate FieldDescr object, for which there is
  /// one per Simulation object (i.e. one per process)

  friend class FieldFace; // required for adjust_alignment_()

public: // interface

  /// Create a new uninitialized FieldData object
  FieldData(const FieldDescr * = 0,
	     int nx=0, int ny=1, int nz=1) throw();

  /// Deconstructor
  ~FieldData() throw();

  /// Copy constructor
  FieldData(const FieldData & field_data) throw ();

  /// Assignment operator
  FieldData & operator= (const FieldData & field_data) throw ();

  void pup(PUP::er &p) ;

  /// Return dimensions of the given field in the block
  void dimensions(const FieldDescr *, int id_field,
		  int * mx, int * my = 0, int * mz = 0) const throw();

  /// Return size of fields on the block, assuming centered
  void size(int * nx, int * ny = 0, int * nz = 0) const throw();

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  char * values (const FieldDescr * , int id_field) 
    throw (std::out_of_range);
  char * values (const FieldDescr * field_descr, std::string name) 
    throw (std::out_of_range)
  { return values (field_descr,field_descr->field_id(name)); }

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  const char * values (const FieldDescr *, int id_field) const 
    throw (std::out_of_range);
  const char * values (const FieldDescr * field_descr, std::string name) const 
    throw (std::out_of_range)
  { return values (field_descr,field_descr->field_id(name)); }

  /// Return array for the corresponding field, which does not contain
  /// ghosts whether they're allocated or not
  char * unknowns (const FieldDescr *, int id_field) 
    throw (std::out_of_range);
  char * unknowns (const FieldDescr * field_descr, std::string name) 
    throw (std::out_of_range)
  { return unknowns (field_descr,field_descr->field_id(name)); }

  const char * unknowns (const FieldDescr *,  int id_field) const 
    throw (std::out_of_range);
  const char * unknowns (const FieldDescr * field_descr, std::string name) const     throw (std::out_of_range)
  { return unknowns (field_descr,field_descr->field_id(name)); }

  /// Return raw pointer to the array of all permanent fields.  Const since
  /// otherwise dangerous due to varying field sizes, precisions,
  /// padding and alignment
  const char * permanent ()  const throw () 
  { return permanent_allocated() ? &array_permanent_[0] : NULL; };

  /// Return width of cells along each dimension
  void cell_width(double xm,   double xp,   double * hx,
		  double ym=0, double yp=0, double * hy=0,
		  double zm=0, double zp=0, double * hz=0) const throw ();

  /// Clear specified array(s) to specified value
  void clear ( const FieldDescr *,
	       float value = 0.0, 
	       int id_field_first = -1, 
	       int id_field_last  = -1) throw();
 
  /// Return whether array is allocated or not
  bool permanent_allocated() const throw()
  { return array_permanent_.size() > 0; }

  /// Return whether array is allocated or not
  size_t permanent_size() const throw()
  { return array_permanent_.size(); }

  /// Allocate storage for the permanent fields
  void allocate_permanent(const FieldDescr *,
			  bool ghosts_allocated = false) throw();

  /// Allocate storage for the temporary fields
  void allocate_temporary(const FieldDescr *,
			  int id) throw (std::out_of_range);

  /// Reallocate storage for the field data, e.g. when changing
  /// from ghosts to non-ghosts [ costly for large blocks ]
  void reallocate_permanent(const FieldDescr * ,
			    bool ghosts_allocated = false) throw();

  /// Deallocate storage for the permanent fields
  void deallocate_permanent() throw();

  /// Deallocate storage for the temporary fields
  void deallocate_temporary(const FieldDescr *,int id) 
    throw (std::out_of_range);

  /// Return whether ghost cells are allocated or not.  
  bool ghosts_allocated() const throw ()
  {  return ghosts_allocated_; }

  /// Return the number of elements (nx,ny,nz) along each axis, and total
  /// number of bytes n
  int field_size (const FieldDescr *,
		  int id_field, int *nx=0, int *ny=0, int *nz=0) const throw();

  /// Print basic field characteristics for debugging
  void print (const FieldDescr *,
	      const char * message,
	      bool use_file = false) const throw();

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

  /// Move (not copy) array to array_permanent_ and offsets to
  /// offsets_
  void restore_permanent_ 
  (const FieldDescr *,
   const char       * array_from,
   std::vector<int> & offsets_from ) throw (std::out_of_range);

  template <class T>
  void print_
  (const T * field,
   const char * field_name,
   const char * message,
   FILE * fp,
   int ixm,int iym,int izm,
   int ixp,int iyp,int izp,
   int nx, int ny, int nz,
   int gx, int gy ,int gz,
   int nxd,int nyd) const;

private: // attributes

  /// Size of fields, assuming centered
  int size_[3];

  /// Single array of permanent fields
  std::vector<char> array_permanent_;

  /// Array of temporary fields
  std::vector<char *> array_temporary_;

  /// Offsets into values_ of the first element of each field
  std::vector<int> offsets_;

  /// Whether ghost values are allocated or not 
  bool ghosts_allocated_;


};   

#endif /* DATA_FIELD_DATA_HPP */
