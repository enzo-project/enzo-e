// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Comm] Declaration of the CommBlock class

#ifndef COMM_COMMBLOCK_HPP
#define COMM_COMMBLOCK_HPP

#ifdef CELLO_TRACE
#define TRACE_ADAPT(MSG)			\
  index_.print(MSG,-1,2,false,simulation());	\
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

//----------------------------------------------------------------------

class CommBlock : public CBase_CommBlock
{
  /// @class    CommBlock
  /// @ingroup  Comm
  ///
  /// @brief [\ref Comm] Handles parallel communication and
  /// synchronization of mesh Blocks

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

  /// Destructor
  virtual ~CommBlock() throw();

  /// Copy constructor
  CommBlock(const CommBlock & block) throw ()
  /// @param     block  Object being copied
  {  copy_(block); }

  /// Assignment operator
  CommBlock & operator = (const CommBlock & block) throw ()
  /// @param     block  Source object of the assignment
  /// @return    The target assigned object
  {  copy_(block);  return *this; }


  //----------------------------------------------------------------------
  // CHARM
  //----------------------------------------------------------------------

  /// Initialize an empty CommBlock
  CommBlock()  { };

  /// Initialize a migrated CommBlock
  CommBlock (CkMigrateMessage *m) : CBase_CommBlock(m) { };

  /// CHARM pupper
  void pup(PUP::er &p);

  //----------------------------------------------------------------------
  // ACCESS METHODS
  //----------------------------------------------------------------------

//----------------------------------------------------------------------
  
  /// Return the Block associated with this CommBlock
  inline Block * block() throw()
  { return block_; };
  inline const Block * block() const throw()   
  { return block_; };

  /// Return the child Block associated with this CommBlock
  inline Block * child_block() throw()   
  { return child_block_; };
  inline const Block * child_block() const throw()  
  { return child_block_; };

  /// Return the index of the root block containing this block
  inline void index_forest (int * ix, int * iy, int * iz) const throw ()
  { index_.array(ix,iy,iz); }

  /// Return the current cycle number
  int cycle() const throw() 
  { return cycle_; };

  /// Return the current time
  double time() const throw()   
  { return time_; };

  /// Return the level in the Hierarchy
  int level() const throw() 
  {  return index_.level(); };

  /// Return the current timestep
  double dt() const throw() 
  { return dt_; };

  /// Return current cell widths
  void cell_width 
  (double * dx, double * dy = 0, double * dz = 0)
  const throw();

  /// Return the current stopping criteria
  bool stop() const throw() 
  { return stop_; };

  /// Return whether this CommBlock is a leaf in the octree forest
  bool is_leaf() const 
  { return is_leaf_; }

  /// Index of the CommBlock
  const Index & index() const 
  { return index_; }

  const int & face_level (const int if3[3]) const
  { return face_level_[IF3(if3)]; }

  const int & face_level_new (const int if3[3]) const
  { return face_level_new_[IF3(if3)]; }

  const int & child_face_level (const int ic3[3], const int if3[3]) const
  { return child_face_level_[ICF3(ic3,if3)]; }

  const int & child_face_level_new (const int ic3[3], const int if3[3]) const
  { return child_face_level_new_[ICF3(ic3,if3)]; }

  void set_face_level (const int if3[3], int level)
  { face_level_[IF3(if3)] = level; }

  void set_face_level_new (const int if3[3], int level)
  { face_level_new_[IF3(if3)] = level; }

  void set_child_face_level (const int ic3[3], const int if3[3], int level)
  { child_face_level_[ICF3(ic3,if3)] = level;  }

  void set_child_face_level_new (const int ic3[3], const int if3[3], int level)
  { child_face_level_new_[ICF3(ic3,if3)] = level; }

  //----------------------------------------------------------------------
  // GENERAL
  //----------------------------------------------------------------------

  Index neighbor_ (const int if3[3], Index * ind = 0) const;

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
  ( int *ix, int *iy, int *iz,  int *nx, int *ny, int *nz ) const;

  /// Return which block faces lie along a domain boundary
  void is_on_boundary (bool boundary[3][2]) const throw();

  /// Determine whether this CommBlock is a leaf and store the result
  void set_leaf()
  {
    bool value = true;
    for (size_t i = 0; i < children_.size(); i++) {
      // deleted children are replaced with thisProxy
      if (children_.at(i) != index_) value = false;    
    }
    is_leaf_ = value;
  }

  void update_levels_ ()
  {
    face_level_ =       face_level_new_;
    child_face_level_ = child_face_level_new_;
  }

  bool is_child_ (const Index & index) const
  { 
    for (size_t i=0; i<children_.size(); i++) {
      if (children_[i] == index) return true;
    }
    return false;
  }

  /// Initialize child face levels given own face levels
  void initialize_child_face_levels_();

  //--------------------------------------------------
  // COMPUTE
  //--------------------------------------------------

  void p_compute_enter()
  {      compute_enter_(); }
  void q_compute_enter()
  {      compute_enter_(); }
  void r_compute_enter(CkReductionMsg * msg)
  {      compute_enter_(); delete msg; }

  void p_compute_exit()
  {      compute_exit_(); }
  void q_compute_exit()
  {      compute_exit_(); }
  void r_compute_exit(CkReductionMsg * msg)
  {      compute_exit_(); delete msg; }

protected:
  void compute_enter_();
  void compute_begin_();
  void compute_exit_();
public:

  //--------------------------------------------------
  // OUTPUT
  //--------------------------------------------------

  void p_output_enter()
  {      output_enter_(); }
  void q_output_enter()
  {      output_enter_(); }
  void r_output_enter(CkReductionMsg * msg)
  {      output_enter_(); delete msg;}

  void p_output_exit()
  {      output_exit_(); }
  void q_output_exit()
  {      output_exit_(); }
  void r_output_exit(CkReductionMsg * msg)
  {      output_exit_(); delete msg; }

  /// Contribute block data to ith output object in the simulation
  void p_output_write (int index_output);

  /// Contribute block data to the Initial input object
  void p_output_read (int index_input = 0)
  {  INCOMPLETE("CommBlock::p_output_read"); }

protected:
  void output_enter_();
  void output_begin_();
  void output_exit_();
public:

  //--------------------------------------------------
  // ADAPT
  //--------------------------------------------------

  void p_adapt_enter() 
  {      adapt_enter_(); }
  void q_adapt_enter() 
  {      adapt_enter_(); }
  void r_adapt_enter(CkReductionMsg * msg) 
  {      adapt_enter_(); delete msg;}

  void p_adapt_next ()
  {      adapt_next_(); }
  void q_adapt_next ()
  {      adapt_next_(); }
  void r_adapt_next (CkReductionMsg * msg)
  {      adapt_next_(); delete msg; }

  void p_adapt_called() 
  {      adapt_called_(); }
  void q_adapt_called() 
  {      adapt_called_(); }
  void r_adapt_called(CkReductionMsg * msg) 
  {      adapt_called_(); delete msg; }

  void p_adapt_end ()  
  {      adapt_end_(); }
  void q_adapt_end ()  
  {      adapt_end_(); }
  void r_adapt_end (CkReductionMsg * msg)  
  {      adapt_end_(); delete msg;}

  void p_adapt_exit() 
  {      adapt_exit_(); }
  void q_adapt_exit() 
  {      adapt_exit_(); }
  void r_adapt_exit(CkReductionMsg * msg) 
  {      adapt_exit_(); delete msg;}


  /// Parent tells child to delete itself
  void p_adapt_delete();
  void p_adapt_recv_neighbor_level 
  (Index index_debug, int ic3[3], int if3[3], int level_now, int level_new);
  void p_adapt_send_child_data
  (int ic3[3],int na, char * array, int nf, int * child_face_level);

  void adapt_send_neighbors_levels(int level);

protected:
  bool do_adapt_();
  void adapt_enter_();
  void adapt_begin_ ();
  void adapt_next_ ();
  void adapt_end_ ();
  void adapt_exit_();
  void adapt_coarsen_();
  void adapt_refine_();
  void adapt_called_();
  int adapt_compute_desired_level_(int level_maximum);
  void adapt_delete_child_(Index index_child);
public:

  //--------------------------------------------------
  // CONTROL AND SYNCHRONIZATION
  //--------------------------------------------------

  /// Syncronize before continuing with next phase
  void control_sync(int phase);

  /// synchronize with count other chares
  void p_control_sync_count(int phase, int count) 
  {      control_sync_count_(phase,count); }

protected:
  void control_sync_neighbor_(int phase);
  void control_sync_count_(int phase, int count);
  void control_call_phase_ (int phase);
public:

  //--------------------------------------------------
  // REFRESH
  //--------------------------------------------------

  /// Refresh ghost zones and apply boundary conditions
  void p_refresh_enter()  
  {      refresh_enter_(); }
  void q_refresh_enter()  
  {      refresh_enter_(); }
  void r_refresh_enter(CkReductionMsg * msg)  
  {      refresh_enter_(); delete msg; }

  /// Exit the refresh phase after QD
  void p_refresh_exit () 
  {      refresh_exit_(); }
  void q_refresh_exit () 
  {      refresh_exit_(); }
  void r_refresh_exit (CkReductionMsg * msg) 
  {      refresh_exit_(); delete msg;  }

  /// Refresh a FieldFace in same, next-coarser, or next-finer level
  void x_refresh_face
  (int n, char buffer[],  int type_refresh, int if3[3], int ic3[3]) 
  {      refresh_face_(n,buffer,type_refresh,if3,ic3); }

  /// Get restricted data from child when it is deleted
  void x_refresh_child (int n, char buffer[],int ic3[3]);

protected:
  void refresh_enter_();
  void refresh_begin_();
  void refresh_exit_ ();
  void refresh_face_
  (int n, char buffer[],  int type_refresh, int if3[3], int ic3[3]);
  void refresh_face (int type_refresh, Index index, int if3[3], int ic3[3]);
public:

  //--------------------------------------------------
  // STOPPING
  //--------------------------------------------------

  /// Entry method after begin_stopping() to call Simulation::r_stopping()
  void r_stopping_compute_timestep(CkReductionMsg * msg);

  /// Enter the stopping phase
  void p_stopping_enter () 
  {      stopping_enter_(); }
  void q_stopping_enter () 
  {      stopping_enter_(); }
  void r_stopping_enter (CkReductionMsg * msg) 
  {      stopping_enter_(); delete msg;  }

  /// Enter the stopping phase
  void p_stopping_exit () 
  {      stopping_exit_(); }
  void q_stopping_exit () 
  {      stopping_exit_(); }
  void r_stopping_exit (CkReductionMsg * msg) 
  {      stopping_exit_(); delete msg;  }

protected:
  void stopping_enter_();
  void stopping_begin_();
  void stopping_exit_();
public:

  /// Exit the stopping phase to exit
  void p_exit () 
  {      exit_(); }
  void q_exit () 
  {      exit_(); }
  void r_exit (CkReductionMsg * msg) 
  {      exit_(); delete msg;  }

protected:
  void exit_();
public:

  //--------------------------------------------------
  // PERFORMANCE
  //--------------------------------------------------

protected:
  /// Start and stop measuring CommBlock-based performance regions
  void performance_start_
  (int index_region, std::string file="", int line=0);
  void performance_stop_
  (int index_region, std::string file="", int line=0);
  void performance_switch_
  (int index_region, std::string file="", int line=0);
public:

  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  /// Set state
  inline void set_state (int cycle, double time, double dt, bool stop)
  { 
    set_cycle(cycle);
    set_time(time); 
    set_dt(dt); 
    set_stop(stop); 
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

  /// Set CommBlock's stopping criteria
  virtual void set_stop (double stop) throw()
  { stop_  = stop; }

  /// Initialize CommBlock
  virtual void initialize () throw()
  {  }

  /// Return the local simulation object
  Simulation * simulation() const;

  //  int count_neighbors() const;

protected: // functions


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

  void debug_faces_(const char * mesg);

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

  /// Current stopping criteria
  bool stop_;

  //--------------------------------------------------

  // /// Indices of all neighboring CommBlocks
  // std::vector<Index> neighbor_index_;

  /// Index of current initialization routine
  int index_initial_;

  /// MESH REFINEMENT

  /// list of child nodes
  std::vector<Index> children_;

  /// Synchronization counter for ghost refresh
  Sync loop_refresh_;

  /// Synchronization counter for ghost refresh
  Sync sync_coarsen_;

  /// Synchronization counter for p_control_sync
  int  count_sync_[SYNC_SIZE];
  int  max_sync_[SYNC_SIZE];

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

  /// Whether CommBlock is a leaf node during adapt phase (stored not
  /// computed to avoid race condition bug #30
  bool is_leaf_;

#ifdef TRACE_MEMORY
  int trace_mem_;
#endif
};

#endif /* COMM_COMMBLOCK_HPP */

