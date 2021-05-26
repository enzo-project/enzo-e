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
		   const EnzoConfig * enzo_config,
		   int max_initial_level) throw();

  /// Constructor
  EnzoInitialHdf5() throw()
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialHdf5);

  /// CHARM++ migration constructor
  EnzoInitialHdf5(CkMigrateMessage *m)
    : Initial (m),
      level_(0)
  {  }

  /// Destructor
  virtual ~EnzoInitialHdf5() throw()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Initialize a Block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

protected: // functions

  int read_dataset_
  (void ** data, int index, Block * block,
   std::string file_name,
   std::string axis_map,
   std::string dataset,
   int *mx, int *my, int *mz,
   int *nx, int *ny, int *nz,
   int *gx, int *gy, int *gz,
   int n4[4], double h4[4],
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

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  // Only initialize Blocks at this level
  int level_;
  
  std::vector < std::string > field_files_;
  std::vector < std::string > field_datasets_;
  std::vector < std::string > field_coords_;
  std::vector < std::string > field_names_;

  std::vector < std::string > particle_files_;
  std::vector < std::string > particle_datasets_;
  std::vector < std::string > particle_coords_;
  std::vector < std::string > particle_types_;
  std::vector < std::string > particle_attributes_;
};

#endif /* ENZO_ENZO_INITIAL_HDF5_HPP */

