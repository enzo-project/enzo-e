// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#ifndef MESH_FACTORY_HPP
#define MESH_FACTORY_HPP

class Hierarchy;
class IoBlock;
class IoFieldData;
class IoParticleData;

class Factory : public PUP::able 
{
  /// @class    Factory
  /// @ingroup  Mesh 
  /// @brief [\ref Mesh] Abstract class for creating concrete Hierarchy,
  /// Patch, and Block objects

public: // interface

  Factory() throw() : PUP::able()
  { TRACE("Factory::Factory()"); }
 
  /// Destructor (must be present to avoid possible vtable link errors)
  virtual ~Factory() throw() { }

  /// CHARM++ function for determining the specific class in the class hierarchy
  PUPable_decl(Factory);

  /// CHARM++ migration constructor for PUP::able

  Factory (CkMigrateMessage *m) : PUP::able(m) 
  { TRACE("Factory::Factory(CkMigrateMessage*)"); }

  /// CHARM++ Pack / Unpack function
  virtual void pup (PUP::er &p);

  /// Create a new Hierarchy [abstract factory design pattern]
  virtual Hierarchy * create_hierarchy 
  ( int rank, int refinement, int max_level) const throw ();

  /// Create an Input / Output accessor object for Block
  virtual IoBlock * create_io_block ( ) const throw();

  /// Create an Input / Output accessor object for a FieldData
  virtual IoFieldData * create_io_field_data 
  ( const FieldDescr * field_descr ) const throw();

  /// Create an Input / Output accessor object for a ParticleData
  virtual IoParticleData * create_io_particle_data 
  ( const ParticleDescr * particle_descr ) const throw();

  /// Create a new CHARM++ Block array
  virtual CProxy_Block create_block_array
  (
   DataMsg * data_msg,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   int num_field_blocks) const throw();

  /// Create a new coarse blocks under the Block array.  For Multigrid
  ///  solvers.  Arguments are the same as create_block_array(), plus
  ///  minimal level min_level < 0
  virtual void create_subblock_array
  (
   DataMsg * data_msg,
   CProxy_Block * block_array,
   int min_level,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   int num_field_blocks) const throw();

  /// Create a new Block
  virtual Block * create_block
  (
   DataMsg * data_msg,
   CProxy_Block * block_array,
   Index index,
   int nx, int ny, int nz,
   int num_field_data,
   int count_adapt,
   int cycle, double time, double dt,
   int narray, char * array, int refresh_type,
   int num_face_level, int * face_level,
   Simulation * simulation = 0
) const throw();

// NEW CODE: See 161206 notes: implementing data objects bound with
// block_array elements
//  
//   virtual void bound_with_block(const CkArrayID &b)
//   {bound_arrays_.push_back(b);};

// protected:

//   std::vector<CkArrayID> bound_arrays_;
  

};

#endif /* MESH_FACTORY_HPP */

