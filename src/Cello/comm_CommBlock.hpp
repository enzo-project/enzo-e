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
class FieldDescr;
class Hierarchy;
class Simulation;

// NOTE: +2 in IC is so indices can be -1, e.g. for child indicies of
// neighboring blocks
// IC  index for child blocks
#define IC(icx,icy,icz)  (((icx)+2)%2) + 2*( (((icy)+2)%2) + 2*(((icz)+2)%2) )
// NC  number of children
#define NC(rank) (1<<(rank))
// IN  index for neighbors
#define IN(axis,face)  ((face) + 2*(axis))
// NC  number of neighbors
#define NN(rank) (2*(rank))

// rank NN NC NI
// 0    0  1  1
// 1    2  2  2=2*1
// 2    4  4  8=4*2
// 3    6  8  24=6*4

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

  /// create a CommBlock with the given block count, lower extent, block
  /// size, and number of field blocks
  CommBlock
  (
#ifndef CONFIG_USE_CHARM
   Simulation * simulation,
#endif
   Index index,
   int nx, int ny, int nz,
   int level,
   int num_field_blocks,
   int count_adapt,
   bool testing=false
) throw();

  /// Initialize an empty CommBlock
  CommBlock()  { };

#ifdef CONFIG_USE_CHARM

  void pup(PUP::er &p)
{
  TRACEPUP;

  CBase_CommBlock::pup(p);

  p | index_;
  p | cycle_;
  p | time_;
  p | dt_;
  p | index_initial_;
  p | level_;
  p | children_;
  p | neighbors_;
  p | niblings_;
  p | count_coarsen_;
  p | count_adapt_;
  p | loop_refresh_;

}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

  /// Initialize a migrated CommBlock
  CommBlock (CkMigrateMessage *m) 
    : CBase_CommBlock(m) { };

#endif


  /// Index of the CommBlock
  Index index() const { return index_; }

#ifdef CONFIG_USE_CHARM

  /// Initialize CommBlock for the simulation.
  void p_initial();

  // /// Call current Initial::enforce() on the block
  // void p_initial_enforce();

  /// Apply the numerical methods on the block
  void p_compute(int cycle, double time, double dt);

  /// Contribute block data to ith output object in the simulation
  void p_write (int index_output);

  /// Contribute block data to the Initial input object
  void p_read (int index_input = 0)
  {  INCOMPLETE("CommBlock::p_read");  }


  //--------------------------------------------------
  // ADAPT
  //--------------------------------------------------

  /// Entry function after prepare() to call Simulation::p_output()
  void p_output(CkReductionMsg * msg);

  /// Begin the adapt phase of one or more adapt steps
  void p_adapt_begin();

  /// Exit the adapt phase after QD
  void q_adapt_end ();

  /// Begin a single adapt step
  void p_adapt (int count_adapt);
  /// End the adapt step after QD
  void q_adapt ();
  //  void adapt();
  void p_refine();
  void coarsen();
  void p_child_can_coarsen(int ic);
  int determine_adapt();
  /// Update the depth of the given child
  void p_print(std::string message) {  
    index().print(message.c_str());
    TRACE2("%s level %d",message.c_str(),level_);
  }

  //--------------------------------------------------
  // REFRESH
  //--------------------------------------------------

  /// Refresh ghost zones and apply boundary conditions
  void p_refresh_begin();

  /// Exit the refresh phase after QD
  void q_refresh_end ();

  /// update ghost zones with given neighbor in same level
  void refresh_same (Index index,
		     int ix, int iy, int iz,
		     int lgx, int lgy, int lgz);
  /// update ghost zones with given neighbor in coarser level
  void refresh_coarse (Index index);

  /// update ghost zones with given neighbor in same level
  void refresh_fine (Index index);

  /// Refresh a FieldFace in the same level
  void x_refresh_same(int n, char buffer[],int fx, int fy, int fz);

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
		      const Hierarchy * hierarchy,
		      int fx, int fy, int fz,
		      int index_field_set = 0) throw();
#endif

  void set_child(Index index);
  void p_set_neighbor(int v3[3]);
  void p_set_nibling(const int v3[3]);

  void delete_child(Index index);
  void p_delete_neighbor(const int v3[3]);
  void p_delete_nibling(const int v3[3]);

  bool is_child (const Index & index);
  bool is_neighbor (const Index & index);
  bool is_nibling (const Index & index);
  
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
  Block * block() throw() { return block_; };
  const Block * block() const throw() { return block_; };

//----------------------------------------------------------------------

  /// Return the index of the root block containing this block
  void index_forest (int * ibx = 0, int * iby = 0, int * ibz = 0) const throw();

  /// Return the name of the block
  std::string name () const throw()
  {
    std::stringstream convert;
    convert << "block_" << id_();
    return convert.str();
  }

  /// Return the size the CommBlock array
  void size_forest (int * nx, int * ny, int * nz) const throw();

  /// Compute the lower extent of the CommBlock in the domain
  void lower(double * xm, double * ym = 0, double * zm = 0) const throw ();

  /// Compute the upper extent of the CommBlock in the domain
  void upper(double * xp, double * yp = 0, double * zp = 0) const throw ();

  /// Return the index of this CommBlock in global coordinates for its level
  void index_global
  ( int *ix, int *iy, int *iz,  
    int *nx, int *ny, int *nz ) const;

  /// Return the current cycle number
  int cycle() const throw() { return cycle_; };

  /// Return the current time
  double time() const throw() { return time_; };

  /// Return the level in the Hierarchy
  int level() const throw() { return level_; };


  /// Return the current timestep
  double dt() const throw() { return dt_; };
 
  /// Return which block faces lie along a domain boundary
  void is_on_boundary (bool boundary[3][2]) throw();

  /// Return whether this CommBlock is a leaf in the octree forest
  bool is_leaf() const
  {
    for (size_t i = 0; i < children_.size(); i++) {
      // deleted children are replaced with thisProxy
      if (children_.at(i) != index_) return false;    
    }
    return true;
  }


public: // virtual functions

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
    TRACE("CommBlock::initialize()\n");
  }

  /// Return the local simulation object
  Simulation * simulation() const;

protected: // functions

  std::string id_ () const throw ()
  {
    char buffer[27];
    int v3[3];
    index_.values(&v3[0],&v3[1],&v3[2]);
    sprintf (buffer,"%08X-%08X-%08X",
	     v3[0],v3[1],v3[2]);
    return buffer;
  }

  void initialize_  ( int nx, int ny, int nz,    bool testing);

  /// Allocate and copy in attributes from give CommBlock
  void copy_(const CommBlock & block) throw();

  /// Return adapt_coarsen, adapt_refine, or adapt_same given two adapt results
  int reduce_adapt_ (int a1, int a2) const throw();

  /// Return the (lower) indices of the CommBlock in the level, 
  /// and the number of indices

#ifdef CONFIG_USE_CHARM

  /// Apply all initial conditions to this Block
  void apply_initial_() throw();

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

  /// Update boundary conditions
  void update_boundary_ ();

#endif

protected: // attributes

#ifndef CONFIG_USE_CHARM
  /// Pointer to parent simulation object (MPI only)
  Simulation * simulation_;
#endif

  ///
  /// Mesh Block that this CommBlock controls
  Block * block_;

  //--------------------------------------------------

  Index index_;

  //--------------------------------------------------

  /// Current cycle number
  int cycle_;

  /// Current time
  double time_;

  /// Current timestep
  double dt_;

  //--------------------------------------------------
  
  /// Index of current initialization routine
  int index_initial_;

  /// MESH REFINEMENT

  /// Mesh refinement level
  int level_;

  /// list of child nodes
  std::vector<Index> children_;

  /// list of neighbor nodes
  std::vector<Index> neighbors_;

  /// list of nibling nodes
  std::vector<Index> niblings_;

  /// Can coarsen only if all children can coarsen
  int count_coarsen_;

  /// Number of adapt steps in the adapt phase
  int count_adapt_;

#ifdef CONFIG_USE_CHARM
  /// Synchronization counter for ghost refresh
  Sync loop_refresh_;
#endif

};

#endif /* COMM_COMMBLOCK_HPP */

