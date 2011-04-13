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
  /// @brief [\ref Method] Abstract class for creating concrete Mesh,
  /// Patch, and Block objects

public: // interface

  /// Create a new Block  [abstract factory design pattern]
  virtual Block * create_block
  (int ix, int iy, int iz,
   int nx, int ny, int nz,
   double xm, double ym, double zm,
   double hx, double hy, double hz,
   int num_field_blocks = 1) const throw();

};

#endif /* ENZO_ENZO_FACTORY_HPP */

