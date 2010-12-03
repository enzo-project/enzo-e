// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_NODE_H
#define MESH_NODE_H

/// @file     mesh_node.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-09-18
/// @brief    Include file for mesh_node[4|2k|3k|16].hpp files

/// @enum face_enum
/// @brief Faces for 2D nodes
enum face_enum {
  face_R = 0,
  face_U = 1,
  face_L = 2,
  face_D = 3,
  face_count };

/// @enum corner_enum
/// @brief Corners for 2D nodes
enum corner_enum {
  corner_UL = 0,
  corner_DL = 1,
  corner_UR = 2,
  corner_DR = 3 };

#endif /* MESH_NODE_H */
