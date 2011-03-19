// $Id: enzo_EnzoOutputImage.hpp 1896 2010-12-03 23:54:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoOutputImage.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Enzo] Declaration for the EnzoOutputImage component

#ifndef ENZO_ENZO_OUTPUT_IMAGE_HPP
#define ENZO_ENZO_OUTPUT_IMAGE_HPP

class EnzoOutputImage : public Output {

  /// @class    EnzoOutputImage
  /// @ingroup  Enzo
  /// @brief [\ref Enzo] class for writing Enzo image fields to images

public: // functions

  /// Create an uninitialized EnzoOutputImage object
  EnzoOutputImage() throw();

  virtual ~EnzoOutputImage() throw();

public: // virtual functions

  /// Write mesh-related data to disk
  virtual void write 
  ( int index, Mesh * mesh, 
    int cycle, double time,
    bool root_call=true, int ix0=0, int iy0=0, int iz0=0) const throw();

  /// Write a patch-related data to disk; may be called by write (Mesh)
  virtual void write 
  ( int index, Patch * patch, Mesh * mesh,
    int cycle, double time, 
    bool root_call=true, int ix0=0, int iy0=0, int iz0=0) const throw();

  /// Write a block-related to disk; may be called by write (Patch)
  virtual void write 
  ( int index, Block * block, Patch * patch, Mesh * mesh, 
    int cycle, double time, 
    bool root_call=true, int ix0=0, int iy0=0, int iz0=0) const throw();

protected: // attributes

};

#endif /* ENZO_ENZO_OUTPUT_IMAGE_HPP */
