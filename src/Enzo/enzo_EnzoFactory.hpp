// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#ifndef ENZO_ENZO_FACTORY_HPP
#define ENZO_ENZO_FACTORY_HPP

class EnzoFactory : public Factory {

  /// @class    EnzoFactory
  /// @ingroup  Enzo
  /// @brief [\ref Enzo] Abstract class for creating concrete Hierarchy,
  /// Patch, and Block objects

public: // interface

  /// CHARM++ constructor
  EnzoFactory() throw() 
  : Factory()
  { TRACE ("EnzoFactory::EnzoFactory()"); }

  PUPable_decl(EnzoFactory);

  EnzoFactory(CkMigrateMessage *m) : Factory (m) {}

  /// CHARM++ Pack / Unpack function
  virtual void pup (PUP::er &p);

  /// Create the Input / Output accessor object for EnzoBlock
  virtual IoBlock * create_io_block () const throw();

  /// Create a new CHARM++ Block array [abstract factory design pattern]
  virtual CProxy_Block create_block_array
  (int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   int num_field_blocks,
   bool testing=false) const throw();

  /// Create a new coarse blocks under the Block array.  For Multigrid
  ///  solvers.  Arguments are the same as create_block_array(), plus
  ///  minimal level min_level < 0
  /// [abstract factory design pattern] 
  virtual CProxy_Block * create_subblock_array
  (int min_level,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   int num_field_blocks,
   bool testing=false) const throw();

  /// Create a new Block  [abstract factory design pattern]
  virtual Block * create_block
  (
   CProxy_Block * block_array,
   Index index,
   int nx, int ny, int nz,
   int num_field_blocks,
   int count_adapt,
   int cycle, double time, double dt,
   int narray, char * array, int op_array,
   int num_face_level, int * face_level,
   bool testing=false,
   Simulation * simulation = 0) const throw();

};

#endif /* ENZO_ENZO_FACTORY_HPP */

