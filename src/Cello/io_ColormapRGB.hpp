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

  /// CHARM++ PUP::able declaration
  PUPable_decl(ColormapRGB);

  /// CHARM++ migration constructor
  ColormapRGB(CkMigrateMessage *m) {}

  /// Constructor
  ColormapRGB(std::vector<double> rgb) throw()
    : Colormap() { rgb_ = rgb; }

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

  /// Apply the colormap to the supplied float array, returning results in kr[],kg[],kb[]
  virtual void apply (double * kr, double * kg, double * kb,
		      int ndx, int ndy, int ndz,
		      int nx,  int ny,  int nz,
		      float * array);

  /// Apply the colormap to the supplied double array, returning results in kr[],kg[],kb[]
  virtual void apply (double * kr, double * kg, double * kb,
		      int ndx, int ndy, int ndz,
		      int nx,  int ny,  int nz,
		      double * array);
  
protected: // attributes


protected: // attributes

  // NOTE: change pup() function whenever attributes change

  std::vector<double> rgb_;
};

#endif /* IO_COLORMAP_RGB_HPP */

