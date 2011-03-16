// $Id: enzo_EnzoMesh.hpp 2071 2011-03-07 23:23:28Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_MESH_HPP
#define ENZO_ENZO_MESH_HPP

/// @file     enzo_EnzoMesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Mar  9 11:34:56 PST 2011
/// @todo     Consider removing (EnzoMesh : Mesh) and (EnzoPatch : Patch) in favor of template approach for Factory Method design pattern (Design Patterns, GoF p.113 )
/// @brief    [\ref Enzo] Declaration of the EnzoMesh class

class EnzoMesh : public Mesh {

  /// @class    EnzoMesh
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Adaptive mesh refinement hierarchy

public: // functions

  /// Initialize an Mesh object
  EnzoMesh(Factory * factory,
	   GroupProcess * group_process,
	   int nx,  int ny,  int nz,
	   int nbx, int nby, int nbz) throw ()
    : Mesh (factory,group_process,nx,ny,nz,nbx,nby,nbz)
  { }

};

#endif /* ENZO_ENZO_MESH_HPP */

