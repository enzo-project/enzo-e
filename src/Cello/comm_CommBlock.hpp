// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Comm] Declaration of the CommBlock class

// #define DEBUG_USE_MAX

#ifndef COMM_COMMBLOCK_HPP
#define COMM_COMMBLOCK_HPP

#ifdef CELLO_TRACE
#define TRACE_ADAPT(MSG)			\
  index_.print(MSG,-1,2);			\
  TRACE1("this = %p",this);
#else
#define TRACE_ADAPT(MSG)			\
  ; 
#endif

class Block;
class FieldFace;
class Factory;
class GroupProcess;
class FieldDescr;
class Hierarchy;
class Simulation;

#include "parallel.def"

// index for child blocks
#define IC3(ic3)  ( ((ic3[0]+2)%2) + 2*( ((ic3[1]+2)%2) + 2*( ((ic3[2]+2)%2) )))

// number of children
#define NC(rank) (1<<(rank))

// index for neighbors (axis,face)
#define IN(axis,face)  ((face) + 2*(axis))

// index for face (ix,iy,iz)
#define IF3(if3)  ((if3[0]+1) + 3*((if3[1]+1) + 3*(if3[2]+1)))

// index for child ic3[] face if3[]
#define ICF3(ic3,if3)  (IF3(if3) + 27*IC3(ic3))

// number of neighbors
#define NN(rank) (2*(rank))

#include <limits>
// const int face_level_unknown = std::numeric_limits<int>::max();
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

#include "mesh.decl.h"
class CommBlock : public CBase_CommBlock
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
   Index index,
   int nx, int ny, int nz,
   int num_field_blocks,
   int num_adapt_steps,
   int cycle, double time, double dt,
   int narray, char * array, int op_array,
   int num_face_level, int * face_level,
   bool testing=false
) throw();

  /// Initialize an empty CommBlock
  CommBlock()  { };

  /// CHARM pupper
  void pup(PUP::er &p);

//----------------------------------------------------------------------

  /// Initialize a migrated CommBlock
  CommBlock (CkMigrateMessage *m) 
    : CBase_CommBlock(m) { };

  /// Index of the CommBlock
  const Index & index() const { return index_; }

  /// Print a hello message for debugging
  void p_hello(std::string message)
  { 
    index_.print(message.c_str());
  }
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

  void p_adapt_mesh() { adapt_mesh(); }
  void adapt_mesh();

  void notify_neighbors(int level);
  void p_get_neighbor_level (Index index_debug,
			     int ic3[3], int if3[3],
			     int level_now, int level_new);

  /// Child notifies parent that it can (or cannot) coarsen
  void p_child_can_coarsen();
  /// Parent notifies children whether it can or cannot coarsen
  void p_parent_can_coarsen();
  /// Child sends restricted data to parent
  void p_get_child_data(int ic3[3],
			int na, char * array,
			int nf, int * child_face_level);

  void refine();
  /// Parent tells child to delete itself
  void p_delete();

  void q_adapt_next ();
  void q_adapt_end ();

  int determine_adapt();
  /// Update the depth of the given child
  void p_print(std::string message) {  
#ifdef DEBUG_ADAPT
    char buffer[255];
    sprintf (buffer,"%s [%p]",message.c_str(),this);
    index().print(buffer);
    TRACE2("%s level %d",buffer,level());
#endif
  }

  const int & face_level (const int if3[3]) const
  {  return face_level_[IF3(if3)];  }
  const int & face_level_new (const int if3[3]) const
  {  return face_level_new_[IF3(if3)];  }

  const int & child_face_level (const int ic3[3], const int if3[3]) const
  {
    check_child_(ic3,"CommBlock::child_face_level()",__FILE__,__LINE__);
    check_face_ (if3,"CommBlock::child_face_level()",__FILE__,__LINE__);
    return child_face_level_[ICF3(ic3,if3)];  }
  const int & child_face_level_new (const int ic3[3], const int if3[3]) const
  {
    check_child_(ic3,"CommBlock::child_face_level()",__FILE__,__LINE__);
    check_face_ (if3,"CommBlock::child_face_level()",__FILE__,__LINE__);
    return child_face_level_new_[ICF3(ic3,if3)];  }

  void set_face_level (const int if3[3], int level)
  {  face_level_[IF3(if3)] = level;  }
  void set_face_level_new (const int if3[3], int level)
  {
#ifdef DEBUG_USE_MAX
    face_level_new_[IF3(if3)] = std::max(face_level_new_[IF3(if3)],level);  
#else
    face_level_new_[IF3(if3)] = level;  
#endif
  }

  void set_child_face_level (const int ic3[3], const int if3[3], int level)
  {
    check_child_(ic3,"CommBlock::set_child_face_level()",__FILE__,__LINE__);
    check_face_ (if3,"CommBlock::set_child_face_level()",__FILE__,__LINE__);
    child_face_level_[ICF3(ic3,if3)] = level;  
  }
  void set_child_face_level_new (const int ic3[3], const int if3[3], int level)
  {
    TRACE7("set_child_face_level_new child %d %d %d  face %d %d %d  level %d",
	   ic3[0],ic3[1],ic3[2],if3[0],if3[1],if3[2],level);
    check_child_(ic3,"CommBlock::set_child_face_level()",__FILE__,__LINE__);
    check_face_ (if3,"CommBlock::set_child_face_level()",__FILE__,__LINE__);
#ifdef DEBUG_USE_MAX
    child_face_level_new_[ICF3(ic3,if3)] = std::max(child_face_level_new_[ICF3(ic3,if3)],level);  
#else
    child_face_level_new_[ICF3(ic3,if3)] = level;
#endif
  }

  void update_levels_ ()
  {
    face_level_ =       face_level_new_;
    child_face_level_ = child_face_level_new_;
  }

  //--------------------------------------------------
  // REFRESH
  //--------------------------------------------------

  /// Refresh ghost zones and apply boundary conditions
  void p_refresh_begin()
  { refresh_begin(); }
  void refresh_begin();

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

  //==================================================

  void delete_child(Index index);
  bool is_child (const Index & index) const;

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  virtual ~CommBlock() throw();

  //----------------------------------------------------------------------

  /// Copy constructor
  CommBlock(const CommBlock & block) throw ()
  /// @param     block  Object being copied
  {  copy_(block); }

  /// Assignment operator
  CommBlock & operator = (const CommBlock & block) throw ()
  /// @param     block  Source object of the assignment
  /// @return    The target assigned object
  {  copy_(block);  return *this; }

  /// Return the Block associated with this CommBlock
  Block * block() throw() { return block_; };
  const Block * block() const throw() { return block_; };

  /// Return the child Block associated with this CommBlock
  Block * child_block() throw() { return child_block_; };
  const Block * child_block() const throw() { return child_block_; };

//----------------------------------------------------------------------

  /// Return the index of the root block containing this block
  
  inline void index_forest (int * ix, int * iy, int * iz) const throw ()
  { index_.array(ix,iy,iz); }

  // void index_forest (int * ibx, int * iby = 0, int * ibz = 0) const throw();

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
  int level() const throw() 
  {  return index_.level(); };

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

  /// Set state
  inline void set_state (int cycle, double time, double dt)
  { set_cycle(cycle); set_time(time); set_dt(dt); }

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

  Index neighbor_ (const int if3[3], Index * ind = 0) const;

  /// Return the CommBlock's desired refinement level based on
  /// local refinement criteria
  int desired_level_(int level_maximum);

  /// Initialize child face levels given own face levels
  void initialize_child_face_levels_();

  /// Delete child after its data has been received
  void delete_child_(Index index_child);

  /// Determine the number of adapt steps (0, 1 or initial_max_level_)
  void get_num_adapt_steps_();

  /// Start and stop measuring CommBlock-based performance regions
  void start_performance_(int index_region, 
		   std::string file="", int line=0);
  void stop_performance_(int index_region, 
		   std::string file="", int line=0);
  void switch_performance_(int index_region, 
		   std::string file="", int line=0);

  /// Return the child adjacent to the given child in the direction of
  /// the given face
  void facing_child_(int jc3[3], const int ic3[3], const int if3[3]) const;

  /// Return limits of faces of the given child corresponding to the
  /// given face of the parent
  // void loop_limits_faces_ 
  // (int ifm3[3], int ifp3[3], const int if3[3], const int ic3[3]) const;

  /// Return whether the given face of the given child and its parent
  /// intersect
  bool child_is_on_face_(const int ic3[3],const int if3[3]) const;

  /// Return the face of the parent corresponding to the given face
  /// of the given child.  Inverse of loop_limits_faces_
  bool parent_face_(int ipf3[3],const int if3[3], const int ic3[3]) const;

  void check_child_(const int ic3[3], const char * msg, const char * file, int line) const 
  {
    ASSERT5 (msg, "child %d %d %d out of range in file %s line %d",
	     ic3[0],ic3[1],ic3[2],file,line,
	     0 <= ic3[0] && ic3[0] <= 1 &&
	     0 <= ic3[1] && ic3[1] <= 1 &&
	     0 <= ic3[2] && ic3[2] <= 1);
  }

  void check_face_(const int if3[3], const char * msg, const char * file, int line) const 
  {
    ASSERT5 (msg, "face %d %d %d out of range in file %s line %d",
	     if3[0],if3[1],if3[2],file,line,
	     -1 <= if3[0] && if3[0] <= 1 &&
	     -1 <= if3[1] && if3[1] <= 1 &&
	     -1 <= if3[2] && if3[2] <= 1);
  }

  void debug_faces_(const char * mesg)
  {
#ifndef DEBUG_ADAPT
    return;
#endif

    TRACE_ADAPT(mesg);
    int if3[3] = {0};
    int ic3[3] = {0};

    for (ic3[1]=1; ic3[1]>=0; ic3[1]--) {
      for (if3[1]=1; if3[1]>=-1; if3[1]--) {

	index_.print(mesg,-1,2,true);

#ifdef CELLO_DEBUG
	char buffer[40];
	sprintf(buffer,"out.debug.%03d",CkMyPe());
	FILE * fp = fopen (buffer,"a");
#endif


	for (if3[0]=-1; if3[0]<=1; if3[0]++) {
#ifdef CELLO_DEBUG
	  fprintf (fp,(ic3[1]==1) ? "%d " : "  ",face_level(if3));
#endif
	  PARALLEL_PRINTF ((ic3[1]==1) ? "%d " : "  ",face_level(if3));
	}
#ifdef CELLO_DEBUG
	fprintf (fp,"| ");
#endif
	PARALLEL_PRINTF ("| ") ;
	for (if3[0]=-1; if3[0]<=1; if3[0]++) {
#ifdef CELLO_DEBUG
	  fprintf (fp,(ic3[1]==1) ? "%d " : "  ",face_level_new(if3));
#endif
	  PARALLEL_PRINTF ((ic3[1]==1) ? "%d " : "  ",face_level_new(if3));
	}
#ifdef CELLO_DEBUG
	fprintf (fp,"| ");
#endif
	PARALLEL_PRINTF ("| ");
	for (ic3[0]=0; ic3[0]<2; ic3[0]++) {
	  for (if3[0]=-1; if3[0]<=1; if3[0]++) {
	    for (if3[0]=-1; if3[0]<=1; if3[0]++) {
#ifdef CELLO_DEBUG
	      fprintf (fp,"%d ",child_face_level(ic3,if3));
#endif
	      PARALLEL_PRINTF ("%d ",child_face_level(ic3,if3));
	    }
	  }
	}
#ifdef CELLO_DEBUG
	fprintf (fp,"| ");
#endif
	PARALLEL_PRINTF ("| ");
	for (ic3[0]=0; ic3[0]<2; ic3[0]++) {
	  for (if3[0]=-1; if3[0]<=1; if3[0]++) {
	    for (if3[0]=-1; if3[0]<=1; if3[0]++) {
#ifdef CELLO_DEBUG
	      fprintf (fp,"%d ",child_face_level_new(ic3,if3));
#endif
	      PARALLEL_PRINTF ("%d ",child_face_level_new(ic3,if3));
	    }
	  }
	}
#ifdef CELLO_DEBUG
	fprintf (fp,"\n");
#endif
	PARALLEL_PRINTF ("\n");
	fflush(stdout);
#ifdef CELLO_DEBUG
	fclose (fp);
#endif

      }
    }
  }

  std::string id_ () const throw ()
  {
    char buffer[27];
    int v3[3];
    index_.values(v3);
    sprintf (buffer,"%08X-%08X-%08X",
	     v3[0],v3[1],v3[2]);
    return buffer;
  }

  /// Allocate and copy in attributes from give CommBlock
  void copy_(const CommBlock & block) throw();

  /// Return the (lower) indices of the CommBlock in the level, 
  /// and the number of indices

  /// Update face_level_[] for given child in refined Commblock
  void refine_face_level_update_ (Index index_child);

  /// Update face_level_[] for coarsened CommBlock
  void coarsen_face_level_update_ (Index index_child);

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

  void loop_limits_refresh_ 
  (int ifacemin[3], int ifacemax[3]) const throw();

  void loop_limits_nibling_ (int ic3m[3], int ic3p[3], const int if3[3]) const throw();

  FieldFace * load_face_(int * n, char ** a,
			 int if3[3], int ic3[3], bool lg3[3],
			 int op_array);
  void store_face_(int n, char * a,
		   int if3[3], int ic3[3], bool lg3[3],
		   int op_array);
  FieldFace * create_face_(int if3[3], int ic3[3],  bool lg3[3],
			   int op_array);

protected: // attributes

  /// Whether block exists
  /// Mesh Block that this CommBlock controls
  Block * block_;

  /// Block for holding restricted child data 
  Block * child_block_;

  //--------------------------------------------------

  /// Index of this CommBlock in the octree forest
  Index index_;

  /// Desired level for the next cycle
  int level_new_;
  
  //--------------------------------------------------

  /// Current cycle number
  int cycle_;

  /// Current time
  double time_;

  /// Current timestep
  double dt_;

  //--------------------------------------------------

  /// Indices of all neighboring CommBlocks
  std::vector<Index> neighbor_index_;

  /// Index of current initialization routine
  int index_initial_;

  /// MESH REFINEMENT

  /// list of child nodes
  std::vector<Index> children_;

  /// Synchronization counter for ghost refresh
  Sync loop_refresh_;

  /// Synchronization counter for ghost refresh
  Sync sync_coarsen_;

  /// current level of neighbors along each face
  std::vector<int> face_level_;

  /// new level of neighbors along each face
  std::vector<int> face_level_new_;

  /// current level of neighbors accumulated from children that can coarsen
  std::vector<int> child_face_level_;

  /// new level of neighbors accumulated from children that can coarsen
  std::vector<int> child_face_level_new_;

  /// Can coarsen only if all children can coarsen
  int count_coarsen_;

  /// Counter for initial mesh creation
  int level_count_;

  /// Number of adapt steps in the adapt phase
  int adapt_step_;

  /// Current adapt value for the block
  int adapt_;

  /// Phase to call after refresh
  int next_phase_;

  /// whether CommBlock has been coarsened and should be deleted
  bool coarsened_;

  /// Whether CommBlock is marked for deletion
  bool delete_;

};

#endif /* COMM_COMMBLOCK_HPP */

