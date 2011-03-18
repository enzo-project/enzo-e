// $Id: mesh_Factory.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_FACTORY_HPP
#define MESH_FACTORY_HPP

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

class Factory {

  /// @class    Factory
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Abstract class for creating concrete Mesh, Mesh, Patch, and Block objects

public: // interface

  /// Create a new Mesh [abstract factory design pattern]
  virtual Mesh * create_mesh
  (GroupProcess * group_process,
   int nx,  int ny,  int nz,
   int nbx, int nby, int nbz) throw ();

  /// Create a new Patch [abstract factory design pattern]
  virtual Patch * create_patch
  (GroupProcess * group_process,
   int nx,   int ny,  int nz,
   int nbx,  int nby, int nbz) throw();

  /// Create a new Block [abstract factory design pattern]
  virtual Block * create_block
  (Patch * patch,
   FieldDescr * field_descr,
   int ix, int iy, int iz,
   int nx, int ny, int nz,
   int num_field_blocks = 1) throw();

private: // functions


private: // attributes


};

#endif /* MESH_FACTORY_HPP */

