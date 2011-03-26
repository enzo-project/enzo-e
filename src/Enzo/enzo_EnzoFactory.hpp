// $Id: method_EnzoFactory.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_FACTORY_HPP
#define ENZO_ENZO_FACTORY_HPP

/// @file     method_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Method] Declaration of the EnzoFactory class

class EnzoFactory : public Factory {

  /// @class    EnzoFactory
  /// @ingroup  Method
  /// @brief    [\ref Method] Abstract class for creating concrete EnzoMesh, EnzoPatch, and EnzoBlock objects

public: // interface

  /// Create a new Mesh  [abstract factory design pattern]
  virtual Mesh       * create_mesh
  (GroupProcess * group_process,
   int nx,  int ny,  int nz,
   int nbx, int nby, int nbz) throw ();

  /// Create a new Patch  [abstract factory design pattern]
  virtual Patch * create_patch
  (GroupProcess * group_process,
   int nx,   int ny,  int nz,
   int nbx,  int nby, int nbz,
   double xm, double ym, double zm,
   double xp, double yp, double zp) throw();

  /// Create a new Block  [abstract factory design pattern]
  virtual Block * create_block
  (Patch * patch,
   FieldDescr * field_descr,
   int ix, int iy, int iz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double xp, double yp, double zp,
   int num_field_blocks = 1) throw();

};

#endif /* ENZO_ENZO_FACTORY_HPP */

