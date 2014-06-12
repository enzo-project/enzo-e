// See LICENSE_CELLO file for license and copyright information

/// @file     io_ColormapRGB.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

ColormapRGB::ColormapRGB() throw ()
  : Colormap()
{
}

//----------------------------------------------------------------------

ColormapRGB::ColormapRGB(const ColormapRGB & colormap_rgb) throw ()
  : Colormap(colormap_rgb)
/// @param     ColormapRGB  Object being copied
{
}

//----------------------------------------------------------------------

ColormapRGB & ColormapRGB::operator= (const ColormapRGB & colormap_rgb) throw ()
/// @param     ColormapRGB  Source object of the assignment
/// @return    The target assigned object
{
  Colormap::operator= (colormap_rgb);

  return *this;
}

void ColormapRGBpup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
}

ColormapRGB::~ColormapRGB() throw ()
{
}

/// CHARM++ Pack / Unpack function
void ColormapRGB::pup (PUP::er &p)
{
  TRACEPUP;
  Colormap::pup(p);
  // NOTE: change this function whenever attributes change
}

//======================================================================

void ColormapRGB::apply (FieldBlock * field_block) 
{
}

void ColormapRGB::color (char * kr, char * kg, char * kb,
			 int    ix, int    iy, int    iz) 
{
}
