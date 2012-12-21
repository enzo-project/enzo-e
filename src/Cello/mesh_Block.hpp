// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Mesh] Declaration of the Block class

#ifndef MESH_BLOCK_HPP
#define MESH_BLOCK_HPP

class Patch;
class Factory;
class GroupProcess;

#ifdef CONFIG_USE_CHARM
#include "mesh.decl.h"
class Block : public CBase_Block
#else
class Block
#endif
{
  /// @class    Block
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Basic serial block of mesh data

  friend class IoBlock;

public: // interface

  /// Create a single Block in a Patch
  Block
  (
   int ibx, int iby, int ibz,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xmp, double ymp, double zmp,
   double xb,  double yb,  double zb,
#ifdef CONFIG_USE_CHARM
   CProxy_Patch proxy_patch,
#endif
   int patch_id,
   int patch_rank,
   int num_field_blocks
) throw();

#ifdef CONFIG_USE_CHARM

  /// Create a Block chare array in a Patch
  Block
  (
   int nbx, int nby, int nbz,
   int nx,  int ny,  int nz,
   double xmp, double ymp, double zmp,
   double xb,  double yb,  double zb,
   CProxy_Patch proxy_patch,
   int patch_id,
   int patch_rank,
   int num_field_blocks
) throw();

#endif

  /// Initialize an empty Block
  Block() { };

#ifdef CONFIG_USE_CHARM

  /// Initialize a migrated Block
  Block (CkMigrateMessage *m) 
    : CBase_Block(m) { };

#endif

  /// Copy constructor
  Block(const Block & block) throw();

  /// Assignment operator
  Block & operator= (const Block & block) throw();

  /// Destructor
  virtual ~Block() throw();

#ifdef CONFIG_USE_CHARM

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er &p);

#endif

#ifdef CONFIG_USE_CHARM

  //==================================================
  // Charm++ Entry Functions
  //==================================================

  /// Initialize block for the simulation.
  void p_initial();

  // /// Call current Initial::enforce() on the block
  // void p_initial_enforce();

  /// Refresh ghost zones and apply boundary conditions
  void p_refresh() { refresh(); }
  void refresh();

  /// Apply the numerical methods on the block
  void p_compute(int cycle, double time, double dt) { compute(); }
  void compute();

  /// Refresh a FieldFace
  void x_refresh(int n, char buffer[],int fx, int fy, int fz);

  /// Contribute block data to ith output object in the simulation
  void p_write (int index_output);

  /// Contribute block data to the Initial input object
  void p_read (int index_input = 0);

  /// Entry function after prepare() to call Simulation::p_output()
  void p_output(CkReductionMsg * msg);

  //--------------------------------------------------
  // enforce boundary, compute timetstep and stopping, reduce to p_output()
  void prepare();

#else /* ! CONFIG_USE_CHARM */

  //==================================================
  // MPI functions
  //==================================================

  /// Refresh ghost data
  void refresh_ghosts(const FieldDescr * field_descr,
		      const Patch * patch,
		      int fx, int fy, int fz,
		      int index_field_set = 0) throw();
#endif

  //==================================================
  // INLINE FUNCTIONS
  //==================================================

  /// Return the ith Field block
  inline const FieldBlock * field_block (int i=0) const throw();

  /// Return the ith Field block
  inline FieldBlock * field_block (int i=0) throw();

  /// Return domain lower extent
  inline void lower(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return domain upper extent
  inline void upper(double * x, double * y = 0, double * z = 0) const throw ();

  /// Return the position of this Block in the containing Patch 
  inline void index_patch (int * ix, int * iy = 0, int * iz = 0) const throw();

  /// Return the index of this Block in the containing Patch 
  inline int index () const throw();

  /// Return the name of the block within its patch, e.g. "block_3"
  inline std::string name () const throw();

  /// Return the name of the parent patch, e.g. "patch_12"
  inline std::string patch_name () const throw();

  /// Return the size the containing Patch
  inline void patch_size (int * nx, int * ny, int * nz) const throw();

  /// Return the current cycle number
  inline int cycle() const throw() { return cycle_; };

  /// Return the current time
  inline double time() const throw() { return time_; };

  /// Return the current timestep
  inline double dt() const throw() { return dt_; };
 
  /// Return which block faces lie along the given boundaries
  void is_on_boundary (double lower[3], 
		       double upper[3],
		       bool boundary[3][2]) throw();

  //==================================================
  // VIRTUAL FUNCTIONS
  //==================================================

  virtual void allocate (const FieldDescr * field_descr) throw();

  /// Set Block's cycle
  virtual void set_cycle (int cycle) throw()
  { cycle_ = cycle;}

  /// Set Block's time
  virtual void set_time (double time) throw()
  { time_  = time; }

  /// Set Block's timestep
  virtual void set_dt (double dt) throw()
  { dt_  = dt; }

  /// Initialize Block
  virtual void initialize () throw()
  { DEBUG ("DEBUG Block::initialize()\n"); }

  //==================================================
protected: // functions
  //==================================================

  /// Allocate and copy in attributes from give Block
  void copy_(const Block & block) throw();

#ifdef CONFIG_USE_CHARM

  /// Determine which faces require boundary updates or communication
  void determine_boundary_
  (
   bool is_boundary[3][2],
   bool * lxm,
   bool * lxp,
   bool * lym,
   bool * lyp,
   bool * lzm,
   bool * lzp
   );

  /// Call Boundary object on specified faces
  void update_boundary_ 
  (
   bool is_boundary[3][2],
   bool lxm,
   bool lxp,
   bool lym,
   bool lyp,
   bool lzm,
   bool lzp
   );

#endif

  //==================================================
protected: // attributes
  //==================================================

#ifdef CONFIG_USE_CHARM

  /// Counter when refreshing faces
  int count_refresh_face_;

  /// Parent patch
  CProxy_Patch proxy_patch_;

#endif

  /// Number of field blocks (required by CHARM++ PUP::er)
  int num_field_blocks_;

  /// Array of field blocks
  std::vector<FieldBlock *> field_block_;

  /// ID of parent patch
  int patch_id_;

  /// Process of parent patch (for Charm++ IO)
  int patch_rank_;

  /// Index into Patch [redundant with CHARM thisIndex.x .y .z]
  int index_[3];

  /// Size of Patch
  int size_[3];

  /// Lower extent of the box associated with the block [computable from Patch]
  double lower_[3];

  /// Upper extent of the box associated with the block [computable from Patch]
  double upper_[3];

  /// Current cycle number
  int cycle_;

  /// Current time
  double time_;

  /// Current timestep
  double dt_;

};

//======================================================================

const FieldBlock * Block::field_block (int i) const throw()
{ 
  return field_block_.at(i); 
}

//----------------------------------------------------------------------

FieldBlock * Block::field_block (int i) throw()
{ 
  return field_block_.at(i); 
}

//----------------------------------------------------------------------

void Block::lower(double * x , double * y , double * z ) const throw ()
{
  if (x) *x = lower_[0];
  if (y) *y = lower_[1];
  if (z) *z = lower_[2];
}

//----------------------------------------------------------------------

void Block::upper(double * x , double * y , double * z ) const throw ()
{
  if (x) *x = upper_[0];
  if (y) *y = upper_[1];
  if (z) *z = upper_[2];
}

//----------------------------------------------------------------------

void Block::index_patch (int * ix, int * iy, int * iz) const throw ()
{
  if (ix) (*ix) = index_[0]; 
  if (iy) (*iy) = index_[1]; 
  if (iz) (*iz) = index_[2]; 
}

//----------------------------------------------------------------------

int Block::index () const throw ()
{
  return index_[0] + size_[0] * (index_[1] + size_[1] * index_[2]);
}

//----------------------------------------------------------------------

std::string Block::name () const throw()
{
  std::stringstream convert;
  convert << "block_" << index();
  return convert.str();
}

//----------------------------------------------------------------------

std::string Block::patch_name () const throw()
{
  std::stringstream convert;
  convert << "patch_" << patch_id_;
  return convert.str();
}

//----------------------------------------------------------------------

void Block::patch_size (int * nx=0, int * ny=0, int * nz=0) const throw ()
{
  if (nx) (*nx)=size_[0]; 
  if (ny) (*ny)=size_[1]; 
  if (nz) (*nz)=size_[2]; 
}
//----------------------------------------------------------------------


#endif /* MESH_BLOCK_HPP */

