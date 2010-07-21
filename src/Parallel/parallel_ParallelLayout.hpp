// $Id: parallel_ParallelLayout.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_ParallelLayout.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the ParallelLayout class

#ifndef PARALLEL_PARALLEL_LAYOUT_HPP
#define PARALLEL_PARALLEL_LAYOUT_HPP

/// @def index_3_to_1
/// @brief compute index i as a function of indices ix,iy,iz.  Inverse
/// operation of index_3_to_1
#define index_3_to_1(ix,iy,iz,nx,ny,nz)		\
  (ix+nx*(iy+ny*iz))

/// @def index_1_to_3
/// @brief compute indices ix,iy,iz as a function of index i.  Inverse
/// operation of index_3_to_1
#define index_1_to_3(i,ix,iy,iz,nx,ny,nz)	\
  ix = i % nx;					\
  iy = (i / nx) % ny;				\
  iz = i / (nx*ny);

#define index_1_to_x(i,nx,ny,nz) (i%nx)
#define index_1_to_y(i,nx,ny,nz) ((i/nx)%ny)
#define index_1_to_z(i,nx,ny,nz) (i/(nx*ny))

enum axis_type {
  axis_x,
  axis_y,
  axis_z };


class ParallelLayout {

  /// @class    ParallelLayout
  /// @ingroup  Parallel
  /// @brief Specify how a Patch is partitioned into process, thread, and computational blacks

public: // interface

  /// (*) Create a ParallelLayout of the given dimensionality, defaulting to serial
  ParallelLayout() throw();

  /// Set how many process blocks the ParallelLayout will be partitioned into 
  void set_process_blocks(int p0, int p1, int p2) throw();

  /// Set how thread blocks each process block is partitioned into 
  void set_thread_blocks(int t0, int t1, int t2) throw();

  /// Set how many compute blocks each thread block is partitioned into  
  void set_data_blocks(int d0, int d1, int d2) throw();

  /// Return the number of processes in the ParallelLayout 
  int process_count () throw();

  /// Return the number of threads per process in the ParallelLayout 
  int thread_count () throw();

  /// Return the number of data blocks per process in the ParallelLayout  
  int data_blocks_per_process () throw();

  /// Return the number of data blocks per thread in the ParallelLayout  
  int data_blocks_per_thread () throw();

  // Neighbor functions

  /// Return whether the given neighbor is inside the ParallelLayout or outside 
  bool neighbor_is_internal (int ip, int it, int id,
			     axis_type axis, int face);

  /// Return the relative process block in the given direction from
  //  the given (ip,thread,data) block
  int neighbor_process (int ip, int it, int id,
			axis_type axis, int face)  throw();

  /// Return the relative thead block in the given direction from the
  /// given (process,thread,data) block
  int neighbor_thread (int ip, int it, int id,
		       axis_type axis, int face) throw();

  /// Return the bounds associated with the given
  /// (process,thread,block) relative to the Patch block 0 < x,y,z < 1
  void box_extent (int ip, int it, int id,
		   double lower_extent[3],    
		   double upper_extent[3]);

  /// Return the index range assigned to the given
  /// (process,thread,block) given a Patch array size
  void array_indices (int ip, int it, int id,
		      int nx, int ny, int nz,
		      int lower_index[3],
		      int upper_index[3]) throw ();
  // Periodic functions

  /// Set whether each axis is periodic or not 
  void set_periodic (axis_type axis, bool periodic);

  /// Return whether each axis is periodic or not 
  bool is_periodic (axis_type axis);

private: // attributes

  ///  number of distributed memory process blocks
  int np_[3];

  /// number of shared memory thread blocks per process
  int nt_[3];

  /// number of compute blocks per thread
  int nd_[3];

  /// periodic
  bool periodic_[3];

private: // functions

  /// Project block indices ip,it,id along given axis
  int neighbor_project_(int ip, int it, int id, axis_type axis, int face);

  /// Return the absolute block indices for the given
  /// (process,thread,block)
  void block_indices_ (int ip, int it, int id, int block_index[3]) throw ();

};
#endif /* PARALLEL_PARALLEL_LAYOUT_HPP */
