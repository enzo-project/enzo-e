// See LICENSE_CELLO file for license and copyright information

/// @file     io_Colormap.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-06-10
/// @brief    [\ref Io] Declaration of the Colormap class

#ifndef IO_COLOR_MAP_HPP
#define IO_COLOR_MAP_HPP

class Colormap {

  /// @class    Colormap
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for representing and applying a color map
  ///
  /// The Colormap base class is used to represent the interface for
  /// mapping floating-point numerical data to colors.  

public: // interface

  /// Constructor
  Colormap() throw() {};

  /// Copy constructor
  Colormap(const Colormap & color_map) throw() 
  { };

  /// Destructor
  ~Colormap() throw() {};

  /// Assignment operator
  Colormap & operator= (const Colormap & color_map) throw()
  { return *this; }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    p | min_;
    p | max_;
  }

  // ----------------------------------------------------------------------

  void set_limit (double min, double max)
  { 
    min_ = min;
    max_ = max; 
  }

  // ----------------------------------------------------------------------

public: // virtual functions

  /// Pre-compute the color map for the FieldBlock
  virtual void apply (FieldBlock * field_block) = 0;

  /// Return pre-computed color (kr,kg,kb) for index (ix,iy,iz)
  virtual void color (char * kr, char * kg, char * kb,
		      int    ix, int    iy, int    iz) = 0;

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  double min_, max_;

};

#endif /* IO_COLOR_MAP_HPP */

