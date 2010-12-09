// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Layout.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the Layout class

#ifndef PARALLEL_LAYOUT_HPP
#define PARALLEL_LAYOUT_HPP

/// @def index_3_to_1
/// @brief compute index i as a function of indices ix,iy,iz.  Inverse
/// operation of index_3_to_1
// #define index_3_to_1(ix,iy,iz,nx,ny,nz)		\
//   (ix+nx*(iy+ny*iz))

/// @def index_1_to_3
/// @brief compute indices ix,iy,iz as a function of index i.  Inverse
/// operation of index_3_to_1
// #define index_1_to_3(i,ix,iy,iz,nx,ny,nz)	\
//   ix = i % nx;					\
//   iy = (i / nx) % ny;				\
//   iz = i / (nx*ny);

// #define index_1_to_x(i,nx,ny,nz) (i%nx)
// #define index_1_to_y(i,nx,ny,nz) ((i/nx)%ny)
// #define index_1_to_z(i,nx,ny,nz) (i/(nx*ny))

enum axis_enum {
  axis_x,
  axis_y,
  axis_z };

enum face_enum {
  face_xm,
  face_xp,
  face_ym,
  face_yp,
  face_zm,
  face_zp };


class Layout {

  /// @class    Layout
  /// @ingroup  Parallel
  /// @brief Specify how a Patch or block is partitioned into blocks or subblocks.

public: // interface

  /// Create a Layout with given "process" start and count
  Layout() throw();

  // initialize / access

  /// Set processor offset and count for the layout
  void set_processors(int process_offset, int process_count) throw();

  /// Set how many blocks in the layout
  void set_blocks(int nb0, int nb1, int nb2) throw();

  /// Return the processor offset and count for the layout
  void processors(int * process_offset, int * process_count) throw();

  /// Return the number of blocks in the layout
  int blocks (int *b0, int *b1, int *b2) throw();

  // Operations

  /// Return the process assigned to the given block
  int process (face_enum face,
	       int ibx, int iby, int ibz)  throw();

@@@@@@@@@@    

private: // attributes

  /// Starting process id
  int process_offset_;

  /// Number of processes
  int process_count_;

  /// number of compute blocks per thread
  int block_count_[3];

private: // functions

};
#endif /* PARALLEL_LAYOUT_HPP */
