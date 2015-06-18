// See LICENSE_CELLO file for license and copyright information

/// @file     io_Colormap.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-06-10
/// @brief    [\ref Io] Declaration of the Colormap class

#ifndef IO_COLOR_MAP_HPP
#define IO_COLOR_MAP_HPP

class Colormap : public PUP::able {

  /// @class    Colormap
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for representing and applying a color map
  ///
  /// The Colormap base class is used to represent the interface for
  /// mapping floating-point numerical data to colors.  

public: // interface

  /// Constructor
  Colormap() throw() 
  : min_(std::numeric_limits<double>::min()),
    max_(std::numeric_limits<double>::max())
  {}

  /// CHARM++ PUP::able declaration
  PUPable_abstract(Colormap);

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    PUP::able::pup(p); 
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

  /// Apply the colormap to the supplied float array, returning results in kr[],kg[],kb[]
  virtual void apply (double * kr, double * kg, double * kb,
		      int ndx, int ndy, int ndz,
		      int nx,  int ny,  int nz,
		      float * array) = 0;

  /// Apply the colormap to the supplied double array, returning results in kr[],kg[],kb[]
  virtual void apply (double * kr, double * kg, double * kb,
		      int ndx, int ndy, int ndz,
		      int nx,  int ny,  int nz,
		      double * array) = 0;

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  double min_, max_;

};

#endif /* IO_COLOR_MAP_HPP */

