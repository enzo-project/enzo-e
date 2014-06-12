// See LICENSE_CELLO file for license and copyright information

/// @file     io_ColormapRGB.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-06-10
/// @brief    [\ref Io] Declaration of the ColormapRGB class

#ifndef IO_COLORMAP_RGB_HPP
#define IO_COLORMAP_RGB_HPP

class FieldBlock;

class ColormapRGB : public Colormap {

  /// @class    ColormapRGB
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for representing an RGB-based color map

public: // interface

  /// Constructor
  ColormapRGB() throw();

  /// Copy constructor
  ColormapRGB(const ColormapRGB & color_mapRGB) throw() ;

  /// Assignment operator
  ColormapRGB & operator= (const ColormapRGB & color_mapRGB) throw();

  /// Destructor
  virtual ~ColormapRGB() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  // ----------------------------------------------------------------------

public: // virtual functions

  /// Pre-compute the color mapRGB for the FieldBlock
  virtual void apply (FieldBlock * field_block);

  /// Return pre-computed color (kr,kg,kb) for index (ix,iy,iz)
  virtual void color (char * kr, char * kg, char * kb,
		      int    ix, int    iy, int    iz);
  
protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Dimensions of the pre-computed color mapping

  int nx,ny,nz;

  /// Red component of the color mapping
  char * r_;
  /// Green component of the color mapping
  char * g_;
  /// Blue component of the color mapping
  char * b_;

};

#endif /* IO_COLORMAP_RGB_HPP */

