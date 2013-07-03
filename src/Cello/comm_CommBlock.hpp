// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Comm] Declaration of the CommBlock class

#ifndef COMM_COMMBLOCK_HPP
#define COMM_COMMBLOCK_HPP

class Block;
class FieldFace;
class Factory;
class GroupProcess;
class FieldDescr;
class Hierarchy;
class Simulation;

// index for child blocks
#define IC(icx,icy,icz)  (((icx)+2)%2) + 2*( (((icy)+2)%2) + 2*(((icz)+2)%2) )

// number of children
#define NC(rank) (1<<(rank))

// index for neighbors (axis,face)
#define IN(axis,face)  ((face) + 2*(axis))

// index for face (ix,iy,iz)
#define IF3(if3)  ((if3[0]+1) + 3*((if3[1]+1) + 3*(if3[2]+1)))

// number of neighbors
#define NN(rank) (2*(rank))

#include <limits>
 const int face_level_unknown = std::numeric_limits<int>::max();
//const int face_level_unknown = -1;

enum phase_type {
  phase_unknown,
  phase_adapt,
  phase_refresh,
  phase_output,
  phase_compute
};

enum array_type {
  op_array_unknown,
  op_array_copy,
  op_array_restrict,
  op_array_prolong
};

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
   int num_field_blocks,
   int count_adapt,
   bool initial,
   int cycle, double time, double dt,
   int narray, char * array, int op_array,
   int num_face_level, int * face_level,
   bool testing=false
) throw();

  /// Initialize an empty CommBlock
  CommBlock()  { };

#ifdef CONFIG_USE_CHARM

  /// CHARM pupper
  void pup(PUP::er &p);

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

  /// Initialize a migrated CommBlock
  CommBlock (CkMigrateMessage *m) 
    : CBase_CommBlock(m) { };

#endif

  /// Index of the CommBlock
  const Index & index() const { return index_; }

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

  /// Begin a single adapt refine step
  void p_adapt_start ();
  /// Start the adapt coarsen step
  void q_adapt_next ();
  /// Stop the adapt step
  void q_adapt_stop ();
  //  void adapt();
  void refine();
  void coarsen();
  void p_balance(int ic3[3], int if3[3],int level);
  bool can_coarsen() const;
  bool can_refine() const;
  void p_child_can_coarsen(int ichild[3],
			   int na, char * array,
			   int nf, int * child_face_level);
  void p_parent_coarsened();
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

  /// send  ghost zones to given neighbor in same level
  void refresh_same (Index index,  int iface[3]);

  /// send ghost zones to coarse neighbor in given direction
  void refresh_coarse (Index index, int iface[3]);

  /// send ghost zones to fine neighbors in given direction
  void refresh_fine (Index index, int iface[3], int ichild[3]);

  /// Refresh a FieldFace in the same level
  void x_refresh_same(int n, char buffer[],  int iface[3]);

  /// Refresh a FieldFace in the next-coarser level
  void x_refresh_coarse(int n, char buffer[],int iface[3],int ichild[3]);

  /// Refresh a FieldFace in the next-finer level
  void x_refresh_fine(int n, char buffer[],  int iface[3],int ichild[3]);

  /// Get restricted data from child when it is deleted
  void x_refresh_child (int n, char buffer[],int ichild[3]);

  //--------------------------------------------------

  /// Output, Monitor, Stopping [reduction], and Timestep [reduction]
  void prepare();

  /// Implementation of refresh
  void refresh();

  //==================================================

#else /* ! CONFIG_USE_CHARM */

  /// Refresh ghost data
  void refresh_ghosts(const FieldDescr * field_descr,
		      const Hierarchy * hierarchy,
		      int fx, int fy, int fz,
		      int index_field_set = 0) throw();
#endif

  void delete_child(Index index);
  bool is_child (const Index & index) const;

  void p_set_face_level (int if3[3], int level, int recurse, int type)
  {
    set_face_level(if3,level,recurse,type); 
  }

  void set_face_level (int if3[3], int level, int recurse, int type);

  const int & face_level (int if3[3]) const
  {  return face_level_[IF3(if3)];  }

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

  /// Return the child Block associated with this CommBlock
  Block * child_block() throw() { return child_block_; };
  const Block * child_block() const throw() { return child_block_; };

//----------------------------------------------------------------------

  /// Return the index of the root block containing this block
  void index_forest (int * ibx, int * iby = 0, int * ibz = 0) const throw();

  /// Return the name of the block
  std::string name () const throw();

  /// Return the size the CommBlock array
  void size_forest (int * nx, int * ny = 0, int * nz = 0) const throw();

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
  void is_on_boundary (bool boundary[3][2]) const throw();

  /// Return whether this CommBlock is a leaf in the octree forest
  bool is_leaf() const
  {
    bool value = true;
    for (size_t i = 0; i < children_.size(); i++) {
      // deleted children are replaced with thisProxy
      if (children_.at(i) != index_) value = false;    
    }
    return value;
  }


public: // virtual functions

  /// Set CommBlock's cycle
  virtual void set_cycle (int cycle) throw()
  { cycle_ = cycle;}

  /// Set CommBlock's time
  virtual void set_time (double time) throw()
  {
    time_  = time; }

  /// Set CommBlock's timestep
  virtual void set_dt (double dt) throw()
  { dt_  = dt; }

  /// Initialize CommBlock
  virtual void initialize () throw()
  {
  }

  /// Return the local simulation object
  Simulation * simulation() const;

protected: // functions

  /// Return the child adjacent to the given child in the direction of
  /// the given face
  void facing_child_(int jc3[3], int ic3[3], int if3[3]) const;

  /// Return limits of faces of the given child corresponding to the
  /// given face of the parent
  void loop_limits_faces_ 
  (int icm3[3], int icp3[3], int if3[3], int ic3[3]) const;

  /// Return whether the given face of the given child and its parent
  /// intersect
  bool child_is_on_face_(int ic3[3],int if3[3]) const;

  /// Return the face of the parent corresponding to the given face
  /// of the given child.  Inverse of loop_limits_faces_
  void parent_face_(int ipf3[3],int if3[3], int ic3[3]) const;

  void debug_faces_(const char *, int *);

  std::string id_ () const throw ()
  {
    char buffer[27];
    int v3[3];
    index_.values(v3);
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

  /// Update face_level_[] for given child in refined Commblock
  void refine_face_level_update_ (Index index_child);

  /// Update face_level_[] for coarsened CommBlock
  void coarsen_face_level_update_ (Index index_child);

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

  void loop_limits_refresh_ 
  (int ifacemin[3], int ifacemax[3]) const throw();

  void loop_limits_nibling_ 
  (int ichildmin[3],
   int ichildmax[3],
   int iface[3]) const throw();

  FieldFace * load_face_(int * narray, char ** array,
			 int iface[3], int ichild[3], bool lghost[3],
			 int op_array);
  FieldFace * store_face_(int narray, char * array,
			 int iface[3], int ichild[3], bool lghost[3],
			 int op_array);
  FieldFace * create_face_(int iface[3], int ichild[3],  bool lghost[3],
			   int op_array);

protected: // attributes

#ifndef CONFIG_USE_CHARM
  /// Pointer to parent simulation object (MPI only)
  Simulation * simulation_;
#endif

  /// Mesh Block that this CommBlock controls
  Block * block_;

  /// Block for holding restricted child data 
  Block * child_block_;

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

#ifdef CONFIG_USE_CHARM
  /// Synchronization counter for ghost refresh
  Sync loop_refresh_;
#endif

  /// level of neighbors (and self) along each face
  std::vector<int> face_level_;

  /// level of neighbors accumulated from children that can coarsen
  std::vector<int> child_face_level_;

  /// Can coarsen only if all children can coarsen
  int count_coarsen_;

  /// Number of adapt steps in the adapt phase
  int count_adapt_;

  /// Current adapt value for the block
  int adapt_;

  /// Phase to call after refresh
  int next_phase_;

  /// whether CommBlock has been coarsened and should be deleted
  bool coarsened_;
  

};

#endif /* COMM_COMMBLOCK_HPP */

