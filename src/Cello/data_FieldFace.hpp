// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldFace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-12
/// @brief    [\ref Data] Interface for the FieldFace class

#ifndef DATA_FIELD_FACE_HPP
#define DATA_FIELD_FACE_HPP

class Prolong;
class Restrict;

class FieldFace {

  /// @class    FieldFace
  /// @ingroup  Data
  /// @brief [\ref Data] Class for loading field faces and storing
  /// field ghosts zones

public: // interface

  static long counter[CONFIG_NODE_SIZE];

  /// Constructor of uninitialized FieldFace

  FieldFace () throw()
  : refresh_type_(refresh_unknown),
    field_list_(),
    prolong_(NULL),
    restrict_(NULL)
  {
    ++counter[cello::index_static()]; 

    for (int i=0; i<3; i++) {
      face_[i] = 0;
      child_[i] = 0;
    }

  }

  /// Constructor of initialized FieldFace
  FieldFace (const Field & field) throw();
     
  /// Destructor
  ~FieldFace() throw();

  /// Copy constructor
  FieldFace(const FieldFace & FieldFace) throw();

  /// Assignment operator
  FieldFace & operator= (const FieldFace & FieldFace) throw();

  /// Comparison operator
  bool operator== (const FieldFace & FieldFace) throw();

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

  /// Set the face
  inline void set_face (int fx, int fy = 0, int fz = 0)
  {
    face_[0] = fx;
    face_[1] = fy;
    face_[2] = fz;
  }
  
  /// Set the face
  inline void face (int * fx, int * fy = 0, int * fz = 0)
  {
    if (fx) (*fx) = face_[0]; 
    if (fy) (*fy) = face_[1]; 
    if (fz) (*fz) = face_[2]; 
  }
  
  inline void invert_face ()
  {
    face_[0] = -face_[0];
    face_[1] = -face_[1];
    face_[2] = -face_[2];
  }
  /// Set child if restrict or prolong is required
  void set_child(int icx, int icy=0, int icz=0) throw()
  {
    child_[0] = icx;
    child_[1] = icy;
    child_[2] = icz;
  }

  /// Set refresh type: refresh_fine(prolong),
  /// refresh_coarse(restrict), or refresh_same(copy)

  void set_refresh (int refresh_type)
  {  refresh_type_ = refresh_type;  }

  /// Set Prolong operation (default is Problem::prolong() )
  void set_prolong (Prolong * prolong)
  { prolong_ = prolong; }
  
  /// Set Restrict operation (default is Problem::restrict() )
  void set_restrict (Restrict * restrict)
  { restrict_ = restrict; }
  
  /// Set the list of fields
  void set_field_list (std::vector<int> const & field_list)
  { field_list_ = field_list; }

  std::vector<int> const & field_list () { return field_list_; }

  /// Create an array with the field's face data
  void face_to_array(Field field, int * n, char ** array) throw();

  /// Use existing array for field's face data
  void face_to_array (Field field, char * array) throw();

  /// Copy the input array data to the field's ghost zones

  void array_to_face(char * array, Field field) throw();

  /// Copy directly the face from source FieldData to destination FieldData
  void face_to_face (Field field_src, Field field_dst);

  /// Calculate the number of bytes needed

  int num_bytes_array (Field field) throw();

  /// Compute loop limits for load or store
  void loop_limits
  (int im3[3], int n3[3], const int nd3[3], const int ng3[3], int refresh_type);

  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer);

  void print (const char * message)
  {
    // char filename[40];
    // sprintf (filename,"ff-%02d.debug",CkMyPe());
    
    //    FILE * fp = fopen (filename,"a");
    //    FILE * fp = stdout;
    CkPrintf (" FieldFace %s %p\n",message,this);
    CkPrintf ("    face_  %d %d %d\n",face_[0],face_[1],face_[2]);
    CkPrintf ("    ghost_  %d %d %d\n",ghost_[0],ghost_[1],ghost_[2]);
    CkPrintf ("    child_  %d %d %d\n",child_[0],child_[1],child_[2]);
    CkPrintf ("    refresh_type %d\n",refresh_type_);
    int n;
    CkPrintf ("    field_list size %d = ",n=field_list_.size());
    for (int i=0; i<n; i++) CkPrintf (" %d",field_list_[i]);
    CkPrintf ("\n");
    //    fclose (fp);
  }
  //--------------------------------------------------

private: // functions

  /// copy data
  void copy_(const FieldFace & field_face); 

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

  /// Precision-agnostic function for copying a field block face
  /// into another block's ghost zones
  template<class T>
  void copy_ (T       * vd, int md3[3], int nd3[3], int id3[3],
	      const T * vs, int ms3[3], int ns3[3], int is3[3]) throw();

private: // attributes

  /// Select face, including edges and corners (-1,-1,-1) to (1,1,1)
  int face_[3];

  /// Whether to include ghost zones along x,y,z axes
  bool ghost_[3];

  /// Child index (0,0,0) to (1,1,1) if restriction or prolongation are used
  int child_[3];

  /// Refresh type: fine, coarse, or same
  int refresh_type_;

  /// List of fields (default all)
  std::vector<int> field_list_;

  /// Prolongation object
  Prolong * prolong_;

  /// Restriction object
  Restrict * restrict_;

};

#endif /* DATA_FIELD_FACE_HPP */
