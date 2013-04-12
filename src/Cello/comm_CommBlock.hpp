// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Comm] Declaration of the CommBlock class

#ifndef COMM_COMMBLOCK_HPP
#define COMM_COMMBLOCK_HPP

class Block;
class Factory;
class GroupProcess;

#ifdef CONFIG_USE_CHARM
#include "mesh.decl.h"
class CommBlock : public CBase_CommBlock
#else
class CommBlock
#endif
{
  /// @class    CommBlock
  /// @ingroup  Comm
  /// @brief    [\ref Comm] Handles parallel communication and synchronization of mesh Blocks

  friend class IoBlock;

public: // interface

  /// create a CommBlock with the given block count, lower forest/patch extent, block
  /// size, and number of field blocks
  CommBlock
  (
   int ibx, int iby, int ibz,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xmp, double ymp, double zmp,
   double xb, double yb, double zb,
#ifdef CONFIG_USE_CHARM
#  ifdef REMOVE_PATCH
#  else
   CProxy_Patch proxy_patch,
#  endif
#endif
   int num_field_blocks
) throw();

#ifdef CONFIG_USE_CHARM

  /// For CHARM CommBlock arrays
  CommBlock
  (
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xmp, double ymp, double zmp,
   double xb, double yb, double zb,
#  ifdef REMOVE_PATCH
#  else
   CProxy_Patch proxy_patch,
#  endif
   int num_field_blocks
) throw();

#endif

  /// Initialize an empty CommBlock
  CommBlock()  { };

#ifdef CONFIG_USE_CHARM

  void pup(PUP::er &p)
{
  TRACEPUP;

  CBase_CommBlock::pup(p);

  p | count_refresh_face_;
#  ifdef REMOVE_PATCH
#  else /* REMOVE_PATCH */
  p | proxy_patch_;
#  endif /* REMOVE_PATCH */

  PUParray(p,index_,3);
  PUParray(p,size_,3);
  PUParray(p,lower_,3);
  PUParray(p,upper_,3);
  p | cycle_;
  p | time_;
  p | dt_;

}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

  /// Initialize a migrated CommBlock
  CommBlock (CkMigrateMessage *m) 
    : CBase_CommBlock(m) { };

#endif

#ifdef CONFIG_USE_CHARM

  /// Initialize CommBlock for the simulation.
  void p_initial();

  // /// Call current Initial::enforce() on the block
  // void p_initial_enforce();

  /// Refresh ghost zones and apply boundary conditions
  void p_refresh() { refresh(); }

  /// Apply the numerical methods on the block
  void p_compute(int cycle, double time, double dt);

  /// Refresh a FieldFace
  void x_refresh(int n, char buffer[],int fx, int fy, int fz);

  /// Contribute block data to ith output object in the simulation
  void p_write (int index_output);

  /// Contribute block data to the Initial input object
  void p_read (int index_input = 0);

  /// Entry function after prepare() to call Simulation::p_output()
  void p_output(CkReductionMsg * msg);

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
#  ifdef REMOVE_PATCH
		      const Hierarchy * hierarchy,
#  else
		      const Patch * patch,
#  endif
		      int fx, int fy, int fz,
		      int index_field_set = 0) throw();
#endif

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  virtual ~CommBlock() throw();

  /// Copy constructor
  CommBlock(const CommBlock & block) throw();

  /// Assignment operator
  CommBlock & operator= (const CommBlock & block) throw();

  /// Return the Block associated with this CommBlock
  Block * block() throw() { return & block_; };
  const Block * block() const throw() { return & block_; };

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

#  ifdef REMOVE_PATCH
  void index_forest (int * ibx = 0, int * iby = 0, int * ibz = 0) const throw();
#  else
  /// Return the position of this CommBlock in the containing Patch 
  void index_patch (int * ibx = 0, int * iby = 0, int * ibz = 0) const throw();
#  endif

  /// Return the index of this CommBlock in the containing Patch / Forest
  int index () const throw();

  /// Return the name of the block
  std::string name () const throw()
  {
    std::stringstream convert;
    convert << "block_" << index();
    return convert.str();
  }

  /// Return the size the containing Patch / Forest
#  ifdef REMOVE_PATCH
  void size_forest (int * nx, int * ny, int * nz) const throw();
#  else
  void size_patch (int * nx, int * ny, int * nz) const throw();
#endif

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

  /// Set CommBlock's cycle
  virtual void set_cycle (int cycle) throw()
  { cycle_ = cycle;}

  /// Set CommBlock's time
  virtual void set_time (double time) throw()
  { time_  = time; }

  /// Set CommBlock's timestep
  virtual void set_dt (double dt) throw()
  { dt_  = dt; }

  /// Initialize CommBlock
  virtual void initialize () throw()
  {
    DEBUG ("DEBUG CommBlock::initialize()\n");
  }

protected: // functions

  void initialize_
  ( int ibx, int iby, int ibz,
    int nbx, int nby, int nbz,
    int nx, int ny, int nz,
    double xpm, double ypm, double zpm, // Domain begin
    double xb, double yb, double zb);    // CommBlock width

  /// Allocate and copy in attributes from give CommBlock
  void copy_(const CommBlock & block) throw();

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

  /// Mesh Block that this CommBlock controls
  Block block_;

#ifdef CONFIG_USE_CHARM

  /// Counter when refreshing faces
  int count_refresh_face_;

#ifdef REMOVE_PATCH

#else
  /// Parent patch
  CProxy_Patch proxy_patch_;
#endif

#endif

  //--------------------------------------------------

  /// Index into Patch / Forest [redundant with CHARM thisIndex.x .y .z]
  int index_[3];

  /// Size of Patch / Forest [redundant with CHARM thisIndex.x .y .z]
  int size_[3];

  /// Lower extent of the box associated with the block [computable]
  double lower_[3];

  /// Upper extent of the box associated with the block [computable]
  double upper_[3];

  //--------------------------------------------------

  /// Current cycle number
  int cycle_;

  /// Current time
  double time_;

  /// Current timestep
  double dt_;

};

#endif /* COMM_COMMBLOCK_HPP */

