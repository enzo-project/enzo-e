// $Id: enzo_EnzoPatch.hpp 2069 2011-03-07 21:27:37Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_PATCH_HPP
#define ENZO_ENZO_PATCH_HPP

/// @file     enzo_EnzoPatch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar  9 11:34:50 PST 2011
/// @brief    [\ref Enzo] Declaration of the interface for the EnzoPatch class

class Factory;

class EnzoPatch : public Patch {

  /// @class    EnzoPatch
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Represent a distributed box of uniform (non-adaptive) data

public: // functions

  /// Constructor for given EnzoPatch size and blocking count
  EnzoPatch(Factory * factory,
	    GroupProcess * group_process,
	    int nx,   int ny,  int nz,
	    int nbx,  int nby, int nbz) throw()
    : Patch (factory,group_process,nx,ny,nz,nbx,nby,nbz)
  { }

};

#endif /* ENZO_ENZO_PATCH_HPP */

