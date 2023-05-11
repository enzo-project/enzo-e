// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Oct 12 14:38:21 PDT 2009
/// @brief    [\ref Data] Fortran-style array class.

#ifndef DATA_FIELD_DATA_HPP
#define DATA_FIELD_DATA_HPP

// #define TRACE_PADDED_FACE

class Block;

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

  friend class Data; // required for set_history_()
  friend class Field;
  friend class FieldFace; // required for adjust_alignment_()

public: // interface

  /// Create a new initialized FieldData object
  FieldData(const FieldDescr * = NULL,
	    int nx=0, int ny=1, int nz=1) throw();

  /// Deconstructor
  ~FieldData() throw();

  void pup(PUP::er &p) ;

  /// Return dimensions of the given field in the block, without assuming that
  /// it is cell-centered. This always includes ghost zones (regardless of
  /// whether they've been allocated).
  void dimensions(const FieldDescr *, int id_field,
		  int * mx, int * my = 0, int * mz = 0) const throw();

  /// Return size of fields on the data, assuming centered (this only includes
  /// the active zone)
  void size(int * nx, int * ny = 0, int * nz = 0) const throw();

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  char * values (const FieldDescr *,
		 int id_field, int history=0) throw ();
  char * values (const FieldDescr * field_descr,
		 std::string name, int history=0) throw ()
  { return values (field_descr,field_descr->field_id(name),history); }

  /// Return array for the corresponding field, which may or may not
  /// contain ghosts depending on if they're allocated
  const char * values (const FieldDescr *,
		       int id_field, int history=0) const throw ();
  const char * values (const FieldDescr * field_descr,
		       std::string name, int history=0) const throw ()
  { return values (field_descr,field_descr->field_id(name),history); }

  /// Return a CelloArray that acts as a view of the corresponding field
  ///
  /// If the field cannot be found the program will abort with an error.
  ///
  /// @tparam T the expected (floating-point) type for the field. If this does
  ///     not match the actual type, the program will abort with an error.
  /// @param field_descr class that stores general field information
  /// @param id_field id specifying the field that is to be loaded
  /// @param choice specifies if ghost zones should be included in the view.
  /// @param history the history index for the specified field
  ///
  /// @returns view of the specified field.
  template<class T>
  CelloArray<T, 3> view(const FieldDescr * field_descr, int id_field,
                        ghost_choice choice = ghost_choice::include,
                        int history=0) throw()
  {
    using noconst_T = typename std::remove_cv<T>::type;
    return make_view_<noconst_T>(field_descr, id_field, choice, history, false);
  }

  template<class T>
  CelloArray<T, 3> view(const FieldDescr * field_descr, std::string name,
                        ghost_choice choice = ghost_choice::include,
                        int history=0) throw()
  { return view<T>(field_descr, field_descr->field_id(name), choice, history); }

  template<class T>
  CelloArray<const T, 3> view(const FieldDescr * field_descr, int id_field,
                              ghost_choice choice = ghost_choice::include,
                              int history=0) const throw()
  {
    return const_cast<FieldData*>(this)->view<T>(field_descr, id_field,
                                                 choice, history);
  }

  template<class T>
  CelloArray<const T, 3> view(const FieldDescr * field_descr, std::string name,
                              ghost_choice choice = ghost_choice::include,
                              int history=0) const throw()
  { return view<T>(field_descr, field_descr->field_id(name), choice, history); }

  /// Return array for the corresponding coarse field
  char * coarse_values (const FieldDescr *,
		 int id_field, int history=0) throw ();
  char * coarse_values (const FieldDescr * field_descr,
		 std::string name, int history=0) throw ()
  { return coarse_values (field_descr,field_descr->field_id(name),history); }

  /// Return array for the corresponding coarse field
  const char * coarse_values (const FieldDescr *,
		       int id_field, int history=0) const throw ();
  const char * coarse_values (const FieldDescr * field_descr,
		       std::string name, int history=0) const throw ()
  { return coarse_values (field_descr,field_descr->field_id(name),history); }

  /// Return a CelloArray that acts as a view of the corresponding coarse field
  ///
  /// If the coarse field cannot be found the program will abort with an error.
  ///
  /// @tparam T the expected (floating-point) type for the field. If this does
  ///     not match the actual type, the program will abort with an error.
  /// @param field_descr class that stores general field information
  /// @param id_field id specifying the field that is to be loaded
  /// @param history the history index for the specified field
  ///
  /// @returns view of the specified coarse field.
  template<class T>
  CelloArray<T, 3> coarse_view(const FieldDescr * field_descr,
                               int id_field, int history=0) throw()
  {
    using noconst_T = typename std::remove_cv<T>::type;
    return make_view_<noconst_T>(field_descr, id_field, ghost_choice::include,
                                 history, true);
  }

  template<class T>
  CelloArray<T, 3> coarse_view(const FieldDescr * field_descr,
                               std::string name, int history=0) throw()
  {
    return coarse_view<T>(field_descr, field_descr->field_id(name), history);
  }

  template<class T>
  CelloArray<const T, 3> coarse_view(const FieldDescr * field_descr,
                                     int id_field, int history=0) const throw()
  {
    return const_cast<FieldData*>(this)->coarse_view<T>(field_descr, id_field,
                                                        history);
  }

  template<class T>
  CelloArray<const T, 3> coarse_view(const FieldDescr * field_descr,
				     std::string name,
                                     int history=0) const throw()
  { return coarse_view<T>(field_descr, field_descr->field_id(name), history); }

  /// Return array for the corresponding field, which does not contain
  /// ghosts whether they're allocated or not
  char * unknowns (const FieldDescr *,
		   int id_field, int history=0) throw ();
  char * unknowns (const FieldDescr * field_descr,
		   std::string name, int history=0) throw ()
  { return unknowns (field_descr,field_descr->field_id(name),history); }

  const char * unknowns (const FieldDescr *,
			 int id_field, int history=0) const throw ();
  const char * unknowns (const FieldDescr * field_descr,
			 std::string name, int history=0) const throw ()
  { return unknowns (field_descr,field_descr->field_id(name),history); }

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
  void allocate_temporary(const FieldDescr *, int id) throw ();

  /// Reallocate storage for the field data, e.g. when changing
  /// from ghosts to non-ghosts [ costly for large blocks ]
  void reallocate_permanent(const FieldDescr * ,
			    bool ghosts_allocated = false) throw();

  /// Deallocate storage for the permanent fields
  void deallocate_permanent() throw();

  /// Deallocate storage for the temporary fields
  void deallocate_temporary(const FieldDescr *,int id) 
    throw ();

  /// Allocate storage for coarse padded array
  void allocate_coarse(const FieldDescr *) throw();
  /// Allocate storage for a coarse padded array
  void allocate_coarse(const FieldDescr *, int id) throw();
  /// Deallocate storage for the coarse fields
  void deallocate_coarse() throw ();
  /// Deallocate storage for the coarse fields
  void deallocate_coarse(int id) throw ();

  /// Return whether ghost cells are allocated or not.  
  bool ghosts_allocated() const throw ()
  {  return ghosts_allocated_; }

  /// Return the number of elements (nx,ny,nz) along each axis
  /// (including ghosts), and total number of bytes n
  int field_size (const FieldDescr *,
		  int id_field, int *nx=0, int *ny=0, int *nz=0) const throw();

  /// Return the number of elements (nx,ny,nz) along each axis of the coarse field
  void coarse_dimensions
  (const FieldDescr *, int id_field, int *nx=0, int *ny=0, int *nz=0) const throw();

  /// Print basic field characteristics for debugging
  void print (const FieldDescr *,
	      const char * message,
	      bool use_file = false) const throw();

  //----------------------------------------------------------------------
  // BLAS Operations [depreciated]
  //----------------------------------------------------------------------

  /// Compute inner product field(ix) . field(iy)
  double dot (const FieldDescr *, int ix, int iy) throw();

  /// Scale vector ix by scalar a
  void scale (const FieldDescr *,
	      int iy, long double a, int ix, bool ghosts = true ) throw();

  //----------------------------------------------------------------------
  // History operations
  //----------------------------------------------------------------------

  /// Copy "current" fields to "old" fields
  void save_history (const FieldDescr *, double time);

  /// Return time for given history
  double history_time (const FieldDescr * field_descr, int ih) const
  {
    const int nh = field_descr->num_history();
    return (1 <= ih && ih <= nh) ? history_time_[ih-1] : 0.0;
  }

  //----------------------------------------------------------------------
  // Units operations
  //----------------------------------------------------------------------

  /// scale the field to cgs units given the unit scaling factor
  /// if it's already in cgs, then leave as-is
  /// except if it's in cgs but the scaling factor has changed (e.g. due to
  /// expansion) then adjust for the new scaling factor
  void units_scale_cgs (const FieldDescr *, int id, double amount);
    
  /// convert the field to "code units" given the unit scaling factor
  /// if it's already in code units, leave it as-is
  /// warning if scaling factor has changed
  void units_scale_code (const FieldDescr *, int id, double amount);

  /// Return the current scaling factor of the given Field
  /// 1.0 if in code units, or the scaling factor if in cgs
  double units_scaling (const FieldDescr *, int id);

  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size (FieldDescr * field_descr) const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (FieldDescr * field_descr, char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (FieldDescr * field_descr, char * buffer);

  //--------------------------------------------------
private: // functions
  //--------------------------------------------------

  template<class T>
  long double dot_(const T* X, const T* Y,
		   int mx, int my, int mz,
		   int nx, int ny, int nz,
		   int gx, int gy, int gz) const throw();

  template<class T>
  void scale_(T* Y, long double a, T * X, bool ghosts,
	      int mx, int my, int mz,
	      int nx, int ny, int nz,
	      int gx, int gy, int gz) const throw();


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
   std::vector<int> & offsets_from ) throw ();

  /// (Re-)initialize temporary fields for history
  void set_history_ (const FieldDescr * field_descr);

  /// Allocate (more) units_scaling_ array values
  void units_allocate_ (int n)
  {
    int i=units_scaling_.size();
    if (i < n+1) {
      units_scaling_.resize(n+1);
      for (; i<n+1; i++) {
	units_scaling_[i] = 1.0;
      }
    }
  }
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


  template<class T>
  CelloArray<T, 3> make_view_
  (const FieldDescr * field_descr,
   int id_field, ghost_choice choice,
   int index_history,  bool coarse) throw();

private: // attributes

  /// Size of fields, assuming centered
  int size_[3];

  /// Single array of permanent fields
  std::vector<char> array_permanent_;

  /// Length of allocated temporary fields
  std::vector<int> temporary_size_;

  /// Array of temporary fields
  std::vector< std::vector<char> > array_temporary_;

  /// Offsets into values_ of the first element of each field
  std::vector<int> offsets_;

  /// Whether ghost values are allocated or not 
  bool ghosts_allocated_;

  /// Temporary field id's used for history.  May be permuted
  /// wrt FieldDescr when copying generations.  Initialized from
  /// FieldDescr copy
  std::vector<int> history_id_;

  /// Saved times for history fields [ip]
  std::vector<double> history_time_;

  /// Current scaling of each field
  std::vector<double> units_scaling_;

  //--------------------------------------------------
  
  /// Length of allocated coarse fields
  std::vector<int> coarse_dimensions_;

  /// Coarse fields with one ghost zone for padded Prolong
  std::vector< std::vector<char> > array_coarse_;

};   

#endif /* DATA_FIELD_DATA_HPP */
