// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#ifndef MESH_FACTORY_HPP
#define MESH_FACTORY_HPP

class GroupProcess;
class Hierarchy;
class IoBlock;
class IoFieldBlock;
class Patch;

class Factory
#ifdef CONFIG_USE_CHARM
  : public PUP::able 
#endif
{
  /// @class    Factory
  /// @ingroup  Mesh 
  /// @brief [\ref Mesh] Abstract class for creating concrete Hierarchy,
  /// Patch, and CommBlock objects

public: // interface

#ifdef CONFIG_USE_CHARM
  Factory() throw() : PUP::able()
  { TRACE("Factory::Factory()"); }
#endif
 
  /// Destructor (must be present to avoid possible vtable link errors)
  virtual ~Factory() throw() { }

#ifdef CONFIG_USE_CHARM

  /// CHARM++ function for determining the specific class in the class hierarchy
  PUPable_decl(Factory);

  /// CHARM++ migration constructor for PUP::able

  Factory (CkMigrateMessage *m) : PUP::able(m) 
  { TRACE("Factory::Factory(CkMigrateMessage*)"); }

  /// CHARM++ Pack / Unpack function
  virtual void pup (PUP::er &p);

#endif

  /// Create a new Hierarchy [abstract factory design pattern]
  virtual Hierarchy * create_hierarchy (int dimension, int refinement) const throw ();

  /// Create a new Patch [abstract factory design pattern]
#ifdef CONFIG_USE_CHARM
  virtual CProxy_Patch * 
#else
  virtual Patch *
#endif
  create_patch 
  (
   const FieldDescr * field_descr,
   int nx,   int ny,  int nz,
   int nx0,  int ny0, int nz0,
   int nbx,  int nby, int nbz,
   double xm, double ym, double zm,
   double xp, double yp, double zp,
   int id,
   bool allocate_blocks = true,
   int process_first=0, int process_last_plus=-1
   ) const throw();

  /// Create an Input / Output accessor object for CommBlock
  virtual IoBlock * create_io_block ( ) const throw();

  /// Create an Input / Output accessor object for a FieldBlock
  virtual IoFieldBlock * create_io_field_block ( ) const throw();

#ifdef CONFIG_USE_CHARM

  /// Create a new CHARM++ CommBlock array
  virtual CProxy_CommBlock create_block_array
  (int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double hx, double hy, double hz,
   CProxy_Patch proxy_patch,
   int patch_id,
   int patch_rank,
   int num_field_blocks = 1,
   bool allocate = true) const throw();

#endif

  /// Create a new CommBlock
  virtual CommBlock * create_block
  (int ibx, int iby, int ibz,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double xb, double yb, double zb,
#ifdef CONFIG_USE_CHARM
   CProxy_Patch proxy_patch,
#endif
   int patch_id,
   int patch_rank,
   int num_field_blocks = 1) const throw();

};

#endif /* MESH_FACTORY_HPP */

