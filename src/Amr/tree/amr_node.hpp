// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef AMR_NODE_H
#define AMR_NODE_H

/// @file     amr_node.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-09-18
/// @brief    Include file for amr_node[4|2k|3k|16].hpp files

/// @enum face_type
/// @brief Faces for 2D nodes: R, U, L, D
enum face_type {
  R = 0,
  U = 1,
  L = 2,
  D = 3,
  num_faces };

/// @enum corner_type
/// @brief Corners for 2D nodes: UL, DL, UR, DR
enum corner_type {
  UL = 0,
  DL = 1,
  UR = 2,
  DR = 3 };

#endif /* AMR_NODE_H */
