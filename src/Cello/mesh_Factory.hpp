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

class Factory {

  /// @class    Factory
  /// @ingroup  Mesh
  /// @brief [\ref Mesh] Abstract class for creating concrete Hierarchy,
  /// Patch, and Block objects

public: // interface

  /// Destructor (must be present to avoid possible vtable link errors)
  virtual ~Factory() throw()
  {}

  /// Create a new Hierarchy [abstract factory design pattern]
  virtual Hierarchy * create_hierarchy () const throw ();

  /// Create a new Patch [abstract factory design pattern]
  virtual Patch * create_patch
  (GroupProcess * group_process,
   int nx,   int ny,  int nz,
   int nx0,  int ny0, int nz0,
   int nbx,  int nby, int nbz,
   double xm, double ym, double zm,
   double xp, double yp, double zp) const throw();

  /// Create an Input / Output accessor object for Block
  virtual IoBlock * create_io_block ( ) const throw();

  /// Create an Input / Output accessor object for a FieldBlock
  virtual IoFieldBlock * create_io_field_block ( ) const throw();

  /// Create a new Block [abstract factory design pattern]
  virtual Block * create_block
  (int ibx, int iby, int ibz,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double xb, double yb, double zb,
   int num_field_blocks = 1) const throw();

#ifdef CONFIG_USE_CHARM
  /// Create a new CHARM++ Block array [abstract factory design pattern]
  virtual CProxy_Block create_block_array
  (int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double hx, double hy, double hz,
   int num_field_blocks = 1) const throw();
#endif

};

#endif /* MESH_FACTORY_HPP */

