// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#ifndef ENZO_ENZO_FACTORY_HPP
#define ENZO_ENZO_FACTORY_HPP

class EnzoMsgCheck;

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

  /// Create a new CHARM++ Block chare array proxy
  virtual CProxy_Block new_block_proxy
  (
   DataMsg * data_msg,
   int nbx, int nby, int nbz) const throw();

  /// Create a new CHARM++ Block array
  virtual void create_block_array
  (
   DataMsg * data_msg,
   CProxy_Block block_array,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   int num_field_blocks) const throw();

  /// Create a new coarse blocks under the Block array.  For Multigrid
  ///  solvers.  Arguments are the same as create_block_array(), plus
  ///  minimal level min_level < 0
  /// [abstract factory design pattern] 
  virtual void create_subblock_array
  (
   DataMsg * data_msg,
   CProxy_Block block_array,
   int min_level,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   int num_field_blocks) const throw();

  /// Create a new Block  [abstract factory design pattern]
  virtual void create_block
  (
   DataMsg * data_msg,
   CProxy_Block block_array,
   Index index,
   int nx, int ny, int nz,
   int num_field_blocks,
   int count_adapt,
   int cycle, double time, double dt,
   int narray, char * array, int refresh_type,
   int num_face_level,
   int * face_level,
   Adapt * adapt,
   Simulation * simulation = 0,
   int io_reader = -1
   ) const throw();

  void create_block_check
  (
   EnzoMsgCheck * msg_check,
   CProxy_Block block_array,
   Index index
   ) const throw();
};

#endif /* ENZO_ENZO_FACTORY_HPP */

