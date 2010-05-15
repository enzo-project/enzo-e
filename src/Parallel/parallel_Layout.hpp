// $Id: parallel_Layout.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Layout.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the Layout class

#ifndef PARALLEL_LAYOUT_HPP
#define PARALLEL_LAYOUT_HPP

#include <vector>

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


class Layout {

  /// @class    Layout
  /// @ingroup  Layout
  /// @brief Specify how a Patch is partitioned into process, thread, and computational blacks

public: // interface

  /// (*) Create a Layout of the given dimensionality, defaulting to serial
  Layout() throw();

  /// ( ) Destructor
  //  ~Layout() throw();

  /// ( ) Copy constructor
  //  Layout(const Layout & layout) throw() throw();

  /// ( ) Assignment operator
  //  Layout & operator= (const Layout & layout) throw() throw();

  /// Set how many process blocks the Layout will be partitioned into 
  void set_process_blocks(int p0, int p1, int p2) throw()
  { np_[0] = p0;
    np_[1] = p1;
    np_[2] = p2; };


  /// Set how thread blocks each process block is partitioned into 
  void set_thread_blocks(int t0, int t1, int t2) throw()
  { nt_[0] = t0;
    nt_[1] = t1;
    nt_[2] = t2; };

  /// Set how many compute blocks each thread block is partitioned into  
  void set_data_blocks(int d0, int d1, int d2) throw()
  { nd_[0] = d0;
    nd_[1] = d1;
    nd_[2] = d2; };

  /// Return the number of processes in the Layout 
  int process_count () throw()
  { return np_[0]*np_[1]*np_[2]; };

  /// Return the number of threads per process in the Layout 
  int thread_count () throw()
  { return  nt_[0]*nt_[1]*nt_[2]; };

  /// Return the number of data blocks per thread in the Layout  
  int data_blocks_per_thread () throw()
  { return nd_[0]*nd_[1]* nd_[2]; }

  /// Return the number of data blocks per process in the Layout  
  int data_blocks_per_process () throw()
  { return thread_count() * data_blocks_per_thread(); };

  // Neighbor functions

  /// Return whether the given neighbor is inside the layout or outside 
  bool neighbor_is_internal (int ip, int it, int id,
			     axis_type axis, int face)
  {

    int k = neighbor_project_(ip,it,id,axis,face);
    return (0 <= k && k < nd_[face] * nt_[face] * np_[face]);
  };

  /// Return the relative process block in the given direction from
  //  the given (ip,thread,data) block
  int neighbor_process (int ip, int it, int id,
			axis_type axis, int face)  throw()
  {

    // Project indices along axis of interest and return updated block index
    int i = neighbor_project_(ip,it,id,axis,face);

    // convert blocks to process blocks, and return relative change
    return index_1_to_z(i,nd_[axis],nt_[axis],np_[axis]) - ip;

  }

  /// Return the relative thead block in the given direction from the
  /// given (process,thread,data) block
  int neighbor_thread (int ip, int it, int id,
		       axis_type axis, int face) throw()
  {


    // Project indices along axis of interest and return updated block index
    int i = neighbor_project_(ip,it,id,axis,face);

    // convert blocks to thread block, and return relative change
    return index_1_to_y(i,nd_[axis],nt_[axis],np_[axis]) - it;

  }

  /// Return the bounds associated with the given
  /// (process,thread,block) relative to the Patch block 0 < x,y,z < 1
  void box_extent (int ip, int it, int id,
		   double lower_extent[3],    
		   double upper_extent[3])
  {

    // Project block id's along each axis
    int ip3[3],it3[3],id3[3];

    index_1_to_3(ip,ip3[0],ip3[1],ip3[2],np_[0],np_[1],np_[2]);
    index_1_to_3(it,it3[0],it3[1],it3[2],nt_[0],nt_[1],nt_[2]);
    index_1_to_3(id,id3[0],id3[1],id3[2],nd_[0],nd_[1],nd_[2]);

    // Convert from (ip,it,id) relative blocks to absolute data blocks
    int kdx,kdy,kdz;
    kdx = index_3_to_1(id3[0],it3[0],ip3[0],nd_[0],nt_[0],np_[0]);
    kdy = index_3_to_1(id3[1],it3[1],ip3[1],nd_[1],nt_[1],np_[1]);
    kdz = index_3_to_1(id3[2],it3[2],ip3[2],nd_[2],nt_[2],np_[2]);

    double hdx = 1.0 / (np_[0]*nt_[0]*nd_[0]);
    double hdy = 1.0 / (np_[1]*nt_[1]*nd_[1]);
    double hdz = 1.0 / (np_[2]*nt_[2]*nd_[2]);

    lower_extent[0] = kdx*hdx;
    lower_extent[1] = kdy*hdy;
    lower_extent[2] = kdz*hdz;

    upper_extent[0] = (kdx+1)*hdx;
    upper_extent[1] = (kdy+1)*hdy;
    upper_extent[2] = (kdz+1)*hdz;

  };

  /// Return the index range assigned to the given
  /// (process,thread,block) given a Patch array size
  void array_indices (int ip, int it, int id,
		      int nx, int ny, int nz,
		      int lower_index[3],
		      int upper_index[3]) throw ()
 
  {
    double lower_extent[3];
    double upper_extent[3];
    box_extent(ip,it,id,lower_extent,upper_extent);

    lower_index[0] = lower_extent[0]*nx;
    lower_index[1] = lower_extent[1]*ny;
    lower_index[2] = lower_extent[2]*nz;
    upper_index[0] = upper_extent[0]*nx;
    upper_index[1] = upper_extent[1]*ny;
    upper_index[2] = upper_extent[2]*nz;

  };

  // Periodic functions

  /// Set whether each axis is periodic or not 
  void set_periodic (axis_type axis, bool periodic)
  { periodic_[axis] = periodic; };

  /// Return whether each axis is periodic or not 
  bool is_periodic (axis_type axis)
  { return periodic_[axis]; };

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

  int neighbor_project_(int ip, int it, int id, axis_type axis, int face)
  {
    int ipa,ita,ida;
    switch (axis) {
    case axis_x:
      ipa = index_1_to_x(ip,np_[0],np_[1],np_[2]);
      ita = index_1_to_x(it,nt_[0],nt_[1],nt_[2]);
      ida = index_1_to_x(id,nd_[0],nd_[1],nd_[2]);
      break;
    case axis_y:
      ipa = index_1_to_y(ip,np_[0],np_[1],np_[2]);
      ita = index_1_to_y(it,nt_[0],nt_[1],nt_[2]);
      ida = index_1_to_y(id,nd_[0],nd_[1],nd_[2]);
      break;
    case axis_z:
      ipa = index_1_to_z(ip,np_[0],np_[1],np_[2]);
      ita = index_1_to_z(it,nt_[0],nt_[1],nt_[2]);
      ida = index_1_to_z(id,nd_[0],nd_[1],nd_[2]);
      break;
    default:
      break;
    }

    // Move "face" blocks (could be < 0) along axis
    ida += face;

    // Convert to total data blocks
    int i = index_3_to_1(ida,ita,ipa,nd_[axis],nt_[axis],np_[axis]);

    // Wrap if periodic
    if (periodic_[axis]) {
      int n = nd_[axis]*nt_[axis]*np_[axis];
      i = (i + n) % n;
    }
    return i;

  };
};
#endif /* PARALLEL_LAYOUT_HPP */
