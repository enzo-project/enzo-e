// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_PATCH_HPP
#define MESH_PATCH_HPP

/// @file     mesh_Patch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @todo     Move "size" to DataBlock's, since that's Field-centric
/// @brief    Declaration of the interface for the Patch class

class Patch {

  /// @class    Patch
  /// @ingroup  Mesh
  /// @brief    Represent a distributed box of uniform (non-adaptive) data

public: // interface

  /// Constructor
  Patch() throw();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Patch() throw();

  /// Copy constructor
  Patch(const Patch & patch) throw();

  /// Assignment operator
  Patch & operator= (const Patch & patch) throw();

  //----------------------------------------------------------------------

  /// Set the data descriptor
  void set_data_descr (DataDescr * data_descr) throw();

  /// Return the data descriptor
  DataDescr * data_descr () const throw();

  /// Set the size of the patch in number of grid cells
  void set_size (int npx, int npy, int npz) throw();

  /// Return the size of the patch in number of grid cells
  void size (int * npx, int * npy=0, int * npz=0) const throw();

  /// Set the layout of the patch, describing processes and blocking
  void set_layout (Layout * layout) throw();

  /// Return the layout of the patch, describing processes and blocking
  Layout * layout () const throw();

  /// Set lower and upper extents of the Patch
  void set_extents (double xm, double xp,
		    double ym, double yp,
		    double zm, double zp) throw();

  /// Return the lower and upper extents of the Patch
  void extents (double * xm,   double * xp,
		double * ym=0, double * yp=0,
		double * zm=0, double * zp=0) const throw();

  
  //--------------------------------------------------

  /// Allocate blocks for given process number
  void allocate(int ip=0) throw();

  /// Deallocate blocks for given process number
  void deallocate(int ip=0) throw();

  /// Whether blocks are allocated for given process number
  bool is_allocated(int ip=0) const throw();

  /// Return the number of data blocks for given process number
  int block_count(int ip=0) const throw();

  /// Return the ith data block for the given process number
  /// WARNING: pointer can only be dereferenced if ip is local
  DataBlock * block(int i, int ip=0) const throw();

  //--------------------------------------------------

public: // entry functions

  /// Initial patch advance, ending with receive_()
  void p_evolve();

  //--------------------------------------------------

private: // attributes

  /// Data descriptor
  DataDescr * data_descr_;

  /// Array of data blocks ib associated with process ip in the patch
  /// data_block[ip][ib]
  typedef std::vector<DataBlock *> data_block_vector;
  std::vector<data_block_vector> data_block_;

  /// Size of the patch
  int size_[3];

  /// Layout: describes blocking, processor range, and block-processor mapping 
  Layout * layout_;

  /// Extent of the patch: xm, xp, ym, yp, zm, zp
  double extents_[6];


};

#endif /* MESH_PATCH_HPP */

