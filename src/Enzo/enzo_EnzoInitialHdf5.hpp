// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialHdf5.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-06-23
/// @brief    [\ref Enzo] Declaration of the EnzoInitialHdf5 class

#ifndef ENZO_ENZO_INITIAL_HDF5_HPP
#define ENZO_ENZO_INITIAL_HDF5_HPP

class EnzoInitialHdf5 : public Initial {

  /// @class    EnzoInitialHdf5
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Read initial conditions from the HDF5 files

public: // interface

  /// Constructor
  EnzoInitialHdf5(int cycle,
		  double time,
		  const EnzoConfig * enzo_config) throw();

  /// Constructor
  EnzoInitialHdf5() throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialHdf5);

  /// CHARM++ migration constructor
  EnzoInitialHdf5(CkMigrateMessage *m)
    : Initial (m)
  {  }

  /// Destructor
  ~EnzoInitialHdf5() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Initialize a Block
  virtual void enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw();

protected: // functions


protected: // attributes

  // NOTE: change pup() function whenever attributes change

  std::vector < std::string > field_files_;
  std::vector < std::string > field_datasets_;
  std::vector < std::string > field_names_;
  std::vector < std::string > field_coords_;

  std::vector < std::string > particle_files_;
  std::vector < std::string > particle_datasets_;
  std::vector < std::string > particle_types_;
  std::vector < std::string > particle_attributes_;

};

#endif /* ENZO_ENZO_INITIAL_HDF5_HPP */

