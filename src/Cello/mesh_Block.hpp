// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Mesh] Declaration of the Block class

#ifndef MESH_BLOCK_HPP
#define MESH_BLOCK_HPP

class FieldDescr;
class FieldBlock;
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

  /// create a Block with the given block count, lower PATCH extent, block
  /// size, and number of field blocks
  Block
  (
   int ibx, int iby, int ibz,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xmp, double ymp, double zmp,
   double xb, double yb, double zb,
#ifdef CONFIG_USE_CHARM
   CProxy_Patch proxy_patch,
#endif
   int patch_id,
   int patch_rank,
   int num_field_blocks
) throw();

#ifdef CONFIG_USE_CHARM

  /// For CHARM Block arrays
  Block
  (
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xmp, double ymp, double zmp,
   double xb, double yb, double zb,
#ifdef CONFIG_USE_CHARM
   CProxy_Patch proxy_patch,
#endif
   int patch_id,
   int patch_rank,
   int num_field_blocks
) throw();

  /// Initialize an empty Block
  Block() { };

  /// Initialize a migrated Block
  Block (CkMigrateMessage *m) { };

#endif

#ifdef CONFIG_USE_CHARM

public: // CHARM++ PUPer

  /// Pack / unpack the Block in a CHARM++ program
  void pup(PUP::er &p);

  /// Initialize block for the simulation.
  void p_initial();

  // /// Call current Initial::enforce() on the block
  // void p_initial_enforce();

  /// Refresh ghost zones and apply boundary conditions
  void p_refresh();

  /// Apply the numerical methods on the block
  void p_compute(int cycle, double time, double dt);

  /// Refresh a FieldFace
  void p_refresh_face(int n, char buffer[],int fx, int fy, int fz);

  /// Contribute block data to ith output object in the simulation
  void p_write (int index_output);

  /// Contribute block data to the Initial input object
  void p_read (int index_input = 0);

  // /// Entry function after initial barrier to call refresh()
  // void p_call_refresh();

  /// Entry function after prepare() to call Simulation::p_output()
  void p_call_output(CkReductionMsg * msg);

  //--------------------------------------------------

  /// Output, Monitor, Stopping [reduction], and Timestep [reduction]
  void prepare();

  /// Implementation of refresh
  void refresh();

  /// Boundary and Method
  void compute();

  //==================================================

#else /* ! CONFIG_USE_CHARM */

  /// Refresh ghost data
  void refresh_ghosts(const FieldDescr * field_descr,
		      const Patch * patch,
		      int fx, int fy, int fz,
		      int index_field_set = 0) throw();
#endif

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  virtual ~Block() throw();

  /// Copy constructor
  Block(const Block & block) throw();

  /// Assignment operator
  Block & operator= (const Block & block) throw();

  //----------------------------------------------------------------------

  /// Return the ith Field block
  const FieldBlock * field_block (int i=0) const throw();

  /// Return the ith Field block
  FieldBlock * field_block (int i=0) throw();

//----------------------------------------------------------------------

  /// Return domain lower extent
  inline void lower(double * x, 
		    double * y = 0,
		    double * z = 0) const throw ()
  {
    if (x) *x = lower_[0];
    if (y) *y = lower_[1];
    if (z) *z = lower_[2];
  }

//----------------------------------------------------------------------

  /// Return domain upper extent
  inline void upper(double * x,
		    double * y = 0,
		    double * z = 0) const throw ()
  {
    if (x) *x = upper_[0];
    if (y) *y = upper_[1];
    if (z) *z = upper_[2];
  }

  /// Return the position of this Block in the containing Patch 
  void index_patch (int * ibx = 0, int * iby = 0, int * ibz = 0) const throw();

  /// Return the index of this Block in the containing Patch 
  int index () const throw();

  /// Return the name of the block within its patch, e.g. "block_3"
  std::string name () const throw()
  {
    std::stringstream convert;
    convert << "block_" << index();
    return convert.str();
  }

  /// Return the name of the parent patch, e.g. "patch_12"
  std::string patch_name () const throw()
  {
    std::stringstream convert;
    convert << "patch_" << patch_id_;
    return convert.str();
  }

  /// Return the size the containing Patch
  void size_patch (int * nx, int * ny, int * nz) const throw();

  /// Return the current cycle number
  int cycle() const throw() { return cycle_; };

  /// Return the current time
  double time() const throw() { return time_; };

  /// Return the current timestep
  double dt() const throw() { return dt_; };
 
  /// Return which block faces lie along the given lower[] and upper[] boundaries
  void is_on_boundary (double lower[3], double upper[3],
		       bool boundary[3][2]) throw();

public: // virtual functions

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
  {
    DEBUG ("DEBUG Block::initialize()\n");
  }

protected: // functions

  /// Allocate and copy in attributes from give Block
  void copy_(const Block & block) throw();

#ifdef CONFIG_USE_CHARM

  /// Determine which faces require boundary updates or communication
  void determine_boundary_
  (
   bool is_boundary[3][2],
   bool * axm,
   bool * axp,
   bool * aym,
   bool * ayp,
   bool * azm,
   bool * azp
   );

  void update_boundary_ 
  (
   bool is_boundary[3][2],
   bool axm,
   bool axp,
   bool aym,
   bool ayp,
   bool azm,
   bool azp
   );

#endif

protected: // attributes

#ifdef CONFIG_USE_CHARM

  /// Counter when refreshing faces
  int count_refresh_face_;

  /// Parent patch
  CProxy_Patch proxy_patch_;

#endif

  //--------------------------------------------------

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

  //--------------------------------------------------

  /// Current cycle number
  int cycle_;

  /// Current time
  double time_;

  /// Current timestep
  double dt_;

};

#endif /* MESH_BLOCK_HPP */

