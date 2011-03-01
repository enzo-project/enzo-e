// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_PATCH_HPP
#define MESH_PATCH_HPP

/// @file     mesh_Patch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @todo     Move "size" to DataBlock's, since that's Field-centric
/// @brief    [\ref Mesh] Declaration of the interface for the Patch class

class Patch {

  /// @class    Patch
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Represent a distributed box of uniform (non-adaptive) data

public: // interface

  /// Constructor for given Patch size and blocking count
  Patch(DataDescr *, 
	int nx,   int ny,  int nz,
	int nbx,  int nby, int nbz) throw();

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

  /// Return the data descriptor
  DataDescr * data_descr () const throw();

  /// Return the size of the patch in number of grid cells
  void size (int * nx, int * ny=0, int * nz=0) const throw();

  /// Return the size of the patch in number of grid cells
  void blocking (int * nbx, int * nby=0, int * nbz=0) const throw();

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

  /// Allocate local blocks
  void allocate_blocks() throw();

  /// Deallocate local blocks
  void deallocate_blocks() throw();

  /// Whether local blocks are allocated
  bool blocks_allocated() const throw();

  /// Return the total number of local blocks
  int num_blocks() const throw();

  /// Return the ith data block
  DataBlock * block(int i) const throw();


  //--------------------------------------------------

public: // entry functions

#ifdef CONFIG_USE_CHARM

  /// Initial patch advance, ending with receive_()
  void p_evolve();

#endif

  //--------------------------------------------------

private: // attributes

  /// Data descriptor
  DataDescr * data_descr_;

  /// Layout: describes blocking, processor range, and block-processor mapping 
  Layout * layout_;

  /// Array of data blocks ib associated with this process
  std::vector<DataBlock * > data_block_;

  /// Size of the patch
  int size_[3];

  /// How the Patch is distributed into DataBlocks
  int blocking_[3];

  /// This process id
  int ip_;

  /// Extent of the patch: xm, xp, ym, yp, zm, zp
  double extents_[6];


};

#endif /* MESH_PATCH_HPP */

