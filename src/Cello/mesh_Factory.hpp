// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#ifndef MESH_FACTORY_HPP
#define MESH_FACTORY_HPP

class Hierarchy;
class IoBlock;
class IoFieldBlock;

class Factory : public PUP::able 
{
  /// @class    Factory
  /// @ingroup  Mesh 
  /// @brief [\ref Mesh] Abstract class for creating concrete Hierarchy,
  /// Patch, and CommBlock objects

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
  (
   int rank, int refinement,
   int process_first, int process_last_plus) const throw ();

  /// Create an Input / Output accessor object for CommBlock
  virtual IoBlock * create_io_block ( ) const throw();

  /// Create an Input / Output accessor object for a FieldBlock
  virtual IoFieldBlock * create_io_field_block ( ) const throw();

  /// Create a new CHARM++ CommBlock array
  virtual CProxy_CommBlock create_block_array
  (
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   int num_field_blocks,
   bool testing = false) const throw();

  /// Create a new CommBlock
  virtual CommBlock * create_block
  (
   CProxy_CommBlock * block_array,
   Index index,
   int nx, int ny, int nz,
   int num_field_blocks,
   int count_adapt,
   int cycle, double time, double dt,
   int narray, char * array, int op_array,
   int num_face_level, int * face_level,
   bool testing = false,
   Simulation * simulation = 0) const throw();

};

#endif /* MESH_FACTORY_HPP */

