// See LICENSE_CELLO file for license and copyright information

/// @file     problem_InitialHdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-06-23
/// @brief    [\ref Problem] Declaration of the InitialHdf5 class

#ifndef PROBLEM_INITIAL_HDF5_HPP
#define PROBLEM_INITIAL_HDF5_HPP

class Config;
class Particle;

class InitialHdf5 : public Initial {

  /// @class    InitialHdf5
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Read initial conditions from the HDF5 HDF5 files

public: // interface

  /// Constructor
  InitialHdf5(int cycle,
		  double time,
                  const Config * config,
                  int max_initial_level) throw();

  /// Constructor
  InitialHdf5() throw()
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(InitialHdf5);

  /// CHARM++ migration constructor
  InitialHdf5(CkMigrateMessage *m)
    : Initial (m),
      level_(0)
  {  }

  /// Destructor
  virtual ~InitialHdf5() throw()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Initialize a Block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

protected: // functions

  template <class T>
  void copy_field_data_to_array_
  (cello_float * array, T * data,
   int mx,int my,int mz,int nx,int ny,int nz,int gx,int gy,int gz,int n4[4],
   int IX, int IY) const;

  template <class T, class S>
  void copy_particle_data_to_array_
  (T * array, S * data,
   Particle particle, int it, int ia, int np);

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

#endif /* PROBLEM_INITIAL_HDF5_HPP */

