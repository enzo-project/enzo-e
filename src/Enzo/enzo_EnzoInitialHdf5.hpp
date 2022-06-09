// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialHdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-05-13
/// @brief    [\ref Enzo] Declaration of the EnzoInitialHdf5 class

#ifndef ENZO_ENZO_INITIAL_HDF5_HPP
#define ENZO_ENZO_INITIAL_HDF5_HPP

class EnzoInitialHdf5 : public Initial {

  /// @class    EnzoInitialHdf5
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Read initial conditions from the HDF5 HDF5 files

public: // interface

  /// Constructor
  EnzoInitialHdf5(int cycle,
		  double time,
                  int max_level,
                  std::string                 format,
                  const int                   blocking[3],
                  std::vector < std::string > field_files,
                  std::vector < std::string > field_datasets,
                  std::vector < std::string > field_coords,
                  std::vector < std::string > field_names,
                  std::vector < std::string > particle_files,
                  std::vector < std::string > particle_datasets,
                  std::vector < std::string > particle_coords,
                  std::vector < std::string > particle_types,
                  std::vector < std::string > particle_attributes
                  ) throw();

  /// Constructor
  EnzoInitialHdf5() throw()
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialHdf5);

  /// CHARM++ migration constructor
  EnzoInitialHdf5(CkMigrateMessage *m)
    : Initial (m),
      max_level_(0),
      i_sync_msg_(-1)
  {  }

  /// Destructor
  virtual ~EnzoInitialHdf5() throw()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initialize a Block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

  void recv_data (Block * block, MsgInitial * msg_initial);

  void copy_dataset_to_field_
  (Block * block,
   std::string field_name, int type_data,
   char * data,
   int mx, int my, int mz,
   int nx, int ny, int nz,
   int gx, int gy, int gz,
   int n4[4], int IX, int IY);

  void copy_dataset_to_particle_
  (Block * block,
   std::string particle_type,
   std::string particle_attribute,
   int type_data,
   char * data,
   int nx, int ny, int nz,
   double h4[4], int IX, int IY, int IZ);

protected: // functions
 
/// Access the Sync counter for messages
  Sync * psync_msg_(Block * block)
  {
    ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
    ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
    return scalar_data->value(scalar_descr,i_sync_msg_);
  }

  char * allocate_array_ (int n, int type_data)
  {
    char * data;
    if (type_data == type_single) {
      data = (char *)new float [n];
    } else if (type_data == type_double) {
      data = (char *)new double [n];
    } else {
      data = nullptr;
      ERROR1 ("EnzoInitialHdf5::allocate_array_()",
              "Unsupported data type %d",type_data);
    }
    return data;
  }
  void delete_array_ (char ** array, int type_data)
  {
    if (type_data == type_single) {
      delete [] (float*)(*array);
    } else if (type_data == type_double) {
      delete [] (double*)(*array);
    } else {
      ERROR1 ("EnzoInitialHdf5::delete_array_()",
              "Unsupported data type %d",type_data);
    }
    (*array) = nullptr;
  }
  int is_reader_ (Index index);
  /// return lower and upper ranges of root blocks overlapping the
  /// given region in the domain
  void root_block_range_ (Index index, int array_lower[3], int array_upper[3]);

  void read_dataset_
  (File * file, char ** data, Index index, int type_data,
   double lower_block[3], double upper_block[3],
   int block_index[3],
   std::string axis_map,
   int nx, int ny, int nz,
   int m4[4], int n4[4], double h4[4],
   int *IX, int *IY, int *IZ);

  template <class T>
  void copy_field_data_to_array_
  (enzo_float * array, T * data,
   int mx,int my,int mz,int nx,int ny,int nz,int gx,int gy,int gz,int n4[4],
   int IX, int IY) const;

  template <class T, class S>
  void copy_particle_data_to_array_
  (T * array, S * data,
   Particle particle, int it, int ia, int np);

  template <class T>
  void update_particle_displacements_
  ( T * array, int nx, int ny, int nz,
    Particle particle, int it, int ia,
    double lower, double h, int axis); 

  void check_cosmology_(File * file) const;
  
protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Only initialize Blocks up to this level
  int max_level_;

  /// Format of the input files: "music", "enzo", "inits"
  typedef std::vector < std::string> vecstr_type;
  
  std::string format_;
  /// Size of the root-level octree array partitioning.  Data in all
  /// blocks within a partition are read from a single root-level
  /// Block
  int         blocking_[3];
  vecstr_type field_files_;
  vecstr_type field_datasets_;
  vecstr_type field_coords_;
  vecstr_type field_names_;

  vecstr_type particle_files_;
  vecstr_type particle_datasets_;
  vecstr_type particle_coords_;
  vecstr_type particle_types_;
  vecstr_type particle_attributes_;
  bool l_particle_displacements_;
  std::string particle_position_names_[3];

  /// Index of the Sync object for counting incoming messages if not
  /// a reader; used to call initial_done() (only) after all expected
  /// messages have been received
  int i_sync_msg_;
};

#endif /* ENZO_ENZO_INITIAL_HDF5_HPP */

