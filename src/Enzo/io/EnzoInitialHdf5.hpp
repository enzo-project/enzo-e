// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialHdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-05-13
/// @brief    [\ref Enzo] Declaration of the EnzoInitialHdf5 class

#ifndef ENZO_ENZO_INITIAL_HDF5_HPP
#define ENZO_ENZO_INITIAL_HDF5_HPP



class DataLoader {

  /// @class    DataLoader
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Helper class used by EnzoInitialHdf5 to
  ///           load initial data from HDF5 files.

  public: // interface

    // Constructor
    DataLoader(Block* block, std::string format);

    // Open the given dataset in the given HDF5 file.
    virtual void open_file(std::string filename, std::string dataset, std::string coordinates);

    // Close any file which may be opened.
    void close_file() {
      file->data_close();
      file->file_close();
    }

    // Load data from an opened file into the block with the given
    // block index.
    void load(int * block_index, Index index_block);

    // Copy data to a local or remote block.
    virtual void copy_data_local(char * data) {}
    virtual void copy_data_remote(Index index_block, char * data) {}

  protected: // methods

    void delete_array_(char ** array, int type_data);
    char * allocate_array_(int n, int type_data);
    void read_dataset_(char ** data, Index index_block, int block_index[3]);
    void check_cosmology_(File * file) const;

    // Attributes
    Block * block;
    FileHdf5 * file;
    int type_data;
    double lower_block[3], upper_block[3];
    int nx, ny, nz, IX, IY, IZ;
    double h4[4];
    int m4[4], n4[4];
    std::string coords, format_;
};



class FieldLoader : public DataLoader {

  /// @class    FieldLoader
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Helper class used to
  ///           load initial field data from HDF5 files.

  public: // interface

    // FieldLoader Constructor
    FieldLoader(Block* block, std::string format) : DataLoader(block, format) {
      Field field = block->data()->field();
      field.ghost_depth(0,&gx,&gy,&gz);
      field.dimensions (0,&mx,&my,&mz);
    }

    // Open the given field dataset in the given HDF5 file.
    void open_file(std::string filename, 
                   std::string field_name,
                   std::string field_coords,
                   std::string field_dataset) {
      DataLoader::open_file(filename, field_dataset, field_coords);
      name = field_name;
    }

    // Load field data from received MsgInitial message.
    void read_msg(MsgInitial * msg_initial, char ** data);
    virtual void copy_data_local(char * data);
    virtual void copy_data_remote(Index index_block, char * data);

  private: // methods

    void copy_dataset_to_field_(char * data);

    template <class T>
    void copy_field_data_to_array_(enzo_float * array, T * data) const;

    // Attributes
    int mx, my, mz, gx, gy, gz;
    std::string name;
};



class ParticleLoader : public DataLoader {

  /// @class    ParticleLoader
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Helper class used to
  ///           load initial particle data from HDF5 files.

  public:

    // ParticleLoader Constructor
    ParticleLoader(Block* block, 
                   bool particle_displacements,
                   std::string format
                   ) : DataLoader(block, format) {
      l_particle_displacements_ = particle_displacements;
    }

    // Open the given particle dataset in the given HDF5 file.
    void open_file(std::string filename,
                  std::string particle_type,
                  std::string particle_attribute,
                  std::string particle_coords,
                  std::string particle_dataset) {
      DataLoader::open_file(filename, particle_dataset, particle_coords);
      type = particle_type;
      attribute = particle_attribute;
    }

    // Load particle data from received MsgInitial message.
    void read_msg(MsgInitial * msg_initial, char ** data);

    virtual void copy_data_local(char * data);
    virtual void copy_data_remote(Index index_block, char * data);

  private: // methods

    void copy_dataset_to_particle_(char * data);

    template <class T, class S>
    void copy_particle_data_to_array_(T * array, S * data, Particle particle, int it, int ia, int np);

    template <class T>
    void update_particle_displacements_(T * array, Particle particle, int it, int ia, double lower, double h, int axis);

    // Attributes
    std::string type, attribute;
    bool l_particle_displacements_;
};



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
                  int                         monitor_iter,
                  std::vector < std::string > field_files,
                  std::vector < std::string > field_datasets,
                  std::vector < std::string > field_coords,
                  std::vector < std::string > field_names,
                  std::vector < int >         field_levels,
                  std::vector < std::string > particle_files,
                  std::vector < std::string > particle_datasets,
                  std::vector < std::string > particle_coords,
                  std::vector < std::string > particle_types,
                  std::vector < std::string > particle_attributes,
                  std::vector < int >         particle_levels
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

  // Get the range of a reader block with index `reader_index` at level `level`.
  void get_reader_range(Index reader_index, int lower[3], int upper[3], int level) throw();

  // Use the provided DataLoader to load data onto blocks at the specified level.
  void load_data(int & count_messages, Block * block, int level, int min_level, DataLoader & loader);

    // initialize particle masses from root level mass constant.
  void initialize_particle_mass(Block * block);

  template <class T>
  void initialize_particle_mass(T * array, Particle particle, int it, int ia, T mass);

protected: // functions
 
/// Access the Sync counter for messages
  Sync * psync_msg_(Block * block)
  {
    ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
    ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
    return scalar_data->value(scalar_descr,i_sync_msg_);
  }

  int is_reader_ (Index index);
  /// return lower and upper ranges of root blocks overlapping the
  /// given region in the domain
  void root_block_range_ (Index index, int array_lower[3], int array_upper[3]);
  
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

  /// Parameter for controling monitoring of progress
  int         monitor_iter_;

  vecstr_type field_files_;
  vecstr_type field_datasets_;
  vecstr_type field_coords_;
  vecstr_type field_names_;
  std::vector < int > field_levels_;

  vecstr_type particle_files_;
  vecstr_type particle_datasets_;
  vecstr_type particle_coords_;
  vecstr_type particle_types_;
  vecstr_type particle_attributes_;
  std::vector < int > particle_levels_;
  bool l_particle_displacements_;
  std::string particle_position_names_[3];

  /// Index of the Sync object for counting incoming messages if not
  /// a reader; used to call initial_done() (only) after all expected
  /// messages have been received
  int i_sync_msg_;
};

#endif /* ENZO_ENZO_INITIAL_HDF5_HPP */

