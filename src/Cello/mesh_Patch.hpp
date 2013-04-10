// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Patch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-05-10
/// @brief    [\ref Mesh] Declaration of the interface for the Patch class

#ifndef MESH_PATCH_HPP
#define MESH_PATCH_HPP

#ifndef CONFIG_USE_CHARM
class Factory;
#endif

#ifdef CONFIG_USE_CHARM
#include "mesh.decl.h"
class Patch : public CBase_Patch
#else
class Patch
#endif
{

  /// @class    Patch
  /// @ingroup  Mesh
  /// @brief [\ref Mesh] Represent a distributed box of uniform
  /// (non-adaptive) data

  friend class IoPatch;

 public: // interface

  /// Constructor for given Patch size and blocking count
  Patch
  (
#ifndef CONFIG_USE_CHARM
   const Factory    * factory,
   const FieldDescr * field_descr,
#endif
   int nx,   int ny,  int nz,
   int nx0,  int ny0, int nz0,
   int nbx,  int nby, int nbz,
   double xm, double ym, double zm,
   double xp, double yp, double zp,
   int id,
   bool allocate_blocks = true,
   int process_first =0, int process_last_plus=-1
   ) throw();


#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// CHARM++ Migration constructor
  Patch(CkMigrateMessage * m) : CBase_Patch(m) 
  { 
    TRACE("Patch::Patch(CkMigrateMessage)");  
  };

#endif

#ifdef CONFIG_USE_CHARM

  /// Wait for all CommBlocks to check in before proceeding
  void s_block(CkCallback function);

  /// Initialize component CommBlocks
  void p_initial ();

  /// Wait for all CommBlocks after p_initial
  void s_initial ();

  /// Wait for all CommBlocks before calling next Output
  void s_output ();

  /// Call write on all CommBlocks
  void p_write(int index_output);

  /// Wait for all CommBlocks after p_write
  void s_write();

  /// Call read on all CommBlocks
  void p_read(int index_input);

  /// Wait for all CommBlocks after p_read
  void s_read();

  /// Apply the numerical methods on the patch
  void p_compute(int cycle, double time, double dt);

  /// Refresh ghost zones and apply boundary conditions
  void p_refresh();

#endif

  /// Destructor
  virtual ~Patch() throw();

  /// Return the size of the patch in number of grid cells
  void size (int * nx, int * ny=0, int * nz=0) const throw();

  /// Return the offset of the patch relative to its parent patch
  void offset (int * nx0, int * ny0=0, int * nz0=0) const throw();

  /// Return the number of CommBlocks along each dimension
  void blocking (int * nbx, int * nby=0, int * nbz=0) const throw();

  /// Return the layout of the patch, describing processes and blocking
  Layout * layout () throw();

  /// Return the layout of the patch, describing processes and blocking
  const Layout * layout () const throw();

  /// Return domain lower extent
  void lower(double * x, double * y=0, double * z=0) const throw ();

  /// Return domain upper extent
  void upper(double * x, double * y=0, double * z=0) const throw ();

  //--------------------------------------------------

  const GroupProcess * group_process()  const throw()
  { return group_process_; };

  //--------------------------------------------------

  /// Return whether CommBlocks have been allocated or not
  bool blocks_allocated() const throw()
  { 
#ifdef CONFIG_USE_CHARM
    return block_exists_;
#else
    return (block_.size() != 0);
#endif
  }

  /// Return the number of CommBlocks
  size_t num_blocks(int * nbx = 0, 
		    int * nby = 0,
		    int * nbz = 0) const throw()
  { 
    if (nbx) *nbx = blocking_[0];
    if (nby) *nby = blocking_[1];
    if (nbz) *nbz = blocking_[2];
    return blocking_[0]*blocking_[1]*blocking_[2]; 
  };

  /// Deallocate local CommBlocks
  void deallocate_blocks() throw();

#ifdef CONFIG_USE_CHARM
  /// Return pointer to the CommBlock CHARM++ chare array
  CProxy_CommBlock * block_array() const throw()
  //  { if (block_exists_) return block_array_; else return 0;}
  { return block_array_;}

#else
    
  /// Return the total number of local blocks
  size_t num_local_blocks() const throw();

  /// Return the ith local CommBlock
  CommBlock * local_block(size_t i) const throw();

#endif

  /// Unit test implementation: to simplify calling from Charm main
  void p_test() throw();

protected: // functions

  /// Allocate array, and optionally allocate element CommBlocks
  void allocate_array_
  (bool allocate_blocks = true,
   const FieldDescr * field_descr = 0) throw ();

protected: // attributes

  /// Array of CommBlocks ib associated with this process
#ifdef CONFIG_USE_CHARM
  CProxy_CommBlock * block_array_;
  bool           block_exists_;
  Loop           block_loop_;
#else
  std::vector<CommBlock * > block_;
#endif

#ifndef CONFIG_USE_CHARM
  /// Factory object for creating CommBlocks
  Factory * factory_;
#endif
  
  /// Parallel Group for distributing the Mesh across processors
  GroupProcess * group_process_;

  /// Layout: describes blocking, processor range, and CommBlock-processor mapping 
  Layout layout_;

  /// Size of the patch
  int size_[3];

  /// Offset of the patch relative to its parent patch 
  int offset_[3];

  /// How the Patch is distributed into CommBlocks
  int blocking_[3];

  /// Lower extent of the patch
  double lower_[3];

  /// Upper extent of the patch
  double upper_[3];

};

#endif /* MESH_PATCH_HPP */

