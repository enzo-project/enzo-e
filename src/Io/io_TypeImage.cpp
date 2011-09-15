// See LICENSE_CELLO file for license and copyright information

/// @file     io_TypeImage.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-12
/// @brief    Implementation of the TypeImage class for writing images

#include "io.hpp"

//----------------------------------------------------------------------

TypeImage::TypeImage() throw ()
{

}

//----------------------------------------------------------------------

void open(File       * file,
	  Simulation * simulation)
{
  // open an existing image type file
  // Allocate and clear image data
  // Initialize image library
}

//----------------------------------------------------------------------



void close()
{
  // write data to disk
  // deallocate and image data
  // finalize image library
}

//----------------------------------------------------------------------



void read(Hierarchy * hierarchy,
	  Patch     * patch,
	  Block     * block,
	  int         field_id)
{
  // Read the field defined on the given block from the file
}

//----------------------------------------------------------------------
  


void write(const Hierarchy * hierarchy,
	   const Patch     * patch,
	   const Block     * block,
	   const int         field_id)
{
  // Write the given field defined on the given block to the file
}

//======================================================================

