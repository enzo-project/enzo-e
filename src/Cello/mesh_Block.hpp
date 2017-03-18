// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Mesh] Declaration of the Block class

#ifndef MESH_BLOCK_HPP
#define MESH_BLOCK_HPP

#ifdef CELLO_TRACE
#define TRACE_ADAPT(MSG)			\
  index_.print(MSG,-1,2,false,simulation());	\
  TRACE1("this = %p",this);
#else
#define TRACE_ADAPT(MSG)			\
  ; 
#endif

class Data;
class MsgRefresh;
class MsgRefine;
class MsgCoarsen;
class Factory;
class FieldDescr;
class FieldFace;
class Hierarchy;
class ItFace;
class ItNeighbor;
class Method;
class Particle;
class ParticleData;
class Refresh;
class Solver;
class Simulation;

//----------------------------------------------------------------------

class Block : public CBase_Block
{
  /// @class    Block
  /// @ingroup  Comm
  ///
  /// @brief [\ref Mesh] Handles parallel communication and
  /// synchronization of mesh Blocks

  friend class IoBlock;

public: // interface

  /// create a Block with the given block count, lower extent, block
  /// size, and number of field blocks
  Block ( MsgRefine * msg );

  // Initialize
  void init (
   Index index,
   int nx, int ny, int nz,
   int num_field_blocks,
   int num_adapt_steps,
   int cycle, double time, double dt,
   int narray, char * array, int refresh_type,
   int num_face_level, int * face_level);

  /// Destructor
  virtual ~Block();

  /// Copy constructor
  Block(const Block & block)
  /// @param     block  Object being copied
  {  copy_(block);  }

  /// Assignment operator
  Block & operator = (const Block & block)
  /// @param     block  Source object of the assignment
  /// @return    The target assigned object
  {  copy_(block);  return *this; }

  //----------------------------------------------------------------------
  // CHARM
  //----------------------------------------------------------------------

  /// Initialize an empty Block
  Block()  { };

  /// Initialize a migrated Block
  Block (CkMigrateMessage *m);

  /// CHARM pupper
  virtual void pup(PUP::er &p);

  //----------------------------------------------------------------------
  // ACCESS METHODS
  //----------------------------------------------------------------------
  
  /// Return the Data associated with this Block
  inline Data * data() throw()
  { return data_; };

  inline const Data * data() const throw()   
  { return data_; };

  /// Return the child Data associated with this Block
  inline Data * child_data() throw()   
  { return child_data_; };
  inline const Data * child_data() const throw()  
  { return child_data_; };

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

  int age() const throw()
  { return age_; };

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

  /// Return whether this Block is a leaf in the octree forest
  bool is_leaf() const 
  { return is_leaf_ && ! (index_.level() < 0); }

  /// Index of the Block
  const Index & index() const 
  { return index_; }

  int face_level (const int if3[3]) const
  { return face_level_curr_[IF3(if3)]; }

  int face_level_next (const int if3[3]) const
  { return face_level_next_[IF3(if3)]; }

  int child_face_level (const int ic3[3], const int if3[3]) const
  { return child_face_level_curr_[ICF3(ic3,if3)]; }

  int child_face_level_next (const int ic3[3], const int if3[3]) const
  { return child_face_level_next_[ICF3(ic3,if3)]; }

  void set_face_level_curr (const int if3[3], int level)
  { face_level_curr_[IF3(if3)] = level; }

  void set_face_level_next (const int if3[3], int level)
  { face_level_next_[IF3(if3)] = level; }

  void set_child_face_level_curr (const int ic3[3], const int if3[3], int level)
  { child_face_level_curr_[ICF3(ic3,if3)] = level;  }

  void set_child_face_level_next (const int ic3[3], const int if3[3], int level)
  { child_face_level_next_[ICF3(ic3,if3)] = level; }

  //----------------------------------------------------------------------
  // GENERAL
  //----------------------------------------------------------------------

  Index neighbor_ (const int if3[3], Index * ind = 0) const;

  /// Return the name of the block
  std::string name () const throw();

  /// Return the size the Block array
  void size_forest (int * nx, int * ny = 0, int * nz = 0) const throw();

  /// Compute the lower extent of the Block in the domain
  void lower(double * xm, double * ym = 0, double * zm = 0) const throw ();

  /// Compute the upper extent of the Block in the domain
  void upper(double * xp, double * yp = 0, double * zp = 0) const throw ();

  /// Return the index of this Block in global coordinates for its level
  void index_global
  ( int *ix, int *iy, int *iz,  int *nx, int *ny, int *nz ) const;

  /// Return which block faces lie along a domain boundary
  void is_on_boundary (bool boundary[3][2]) const throw();

  /// Return which faces are periodic
  void periodicity (bool periodic[3][2]) const;

  void update_levels_ ()
  {
    face_level_curr_ =       face_level_next_;
    //    for (int i=0; i<face_level_next_.size(); i++) face_level_next_[i]=0;
    child_face_level_curr_ = child_face_level_next_;
    //    for (int i=0; i<child_face_level_next_.size(); i++) child_face_level_next_[i]=0;
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

  /// Return an iterator over faces

  ItFace it_face(int min_face_rank,
		 Index index,
		 const int * ic3=0,
		 const int * if3=0) throw();

  /// Return an iterator over neighbors

  ItNeighbor it_neighbor(int min_face_rank, Index index) throw();

  //--------------------------------------------------
  // Charm++ virtual
  //--------------------------------------------------

  virtual const CProxy_Block proxy_array() const 
  { return thisProxy; }

  virtual const CProxyElement_Block proxy_element() const 
  { return thisProxy[thisIndex]; }

  //--------------------------------------------------
  // INITIAL
  //--------------------------------------------------

  void initial_exit_();
  void p_initial_exit()
  {
    performance_start_(perf_initial);
    initial_exit_();
    performance_stop_(perf_initial);
  }
  void r_initial_exit(CkReductionMsg * msg)
  {
    performance_start_(perf_initial);
    initial_exit_(); delete msg;
    performance_stop_(perf_initial);
  }

  //--------------------------------------------------
  // COMPUTE
  //--------------------------------------------------

  void p_compute_enter()
  {
    performance_start_(perf_compute,__FILE__,__LINE__);
    compute_enter_();
    performance_stop_(perf_compute,__FILE__,__LINE__);
  }
  void r_compute_enter(CkReductionMsg * msg)
  {
    performance_start_(perf_compute,__FILE__,__LINE__);
    delete msg;    
    compute_enter_();
    performance_stop_(perf_compute,__FILE__,__LINE__);
  }

  void p_compute_continue()
  {
    performance_start_(perf_compute,__FILE__,__LINE__);
    compute_continue_();
    performance_stop_(perf_compute,__FILE__,__LINE__);
  }
  void r_compute_continue(CkReductionMsg * msg)
  {
    performance_start_(perf_compute,__FILE__,__LINE__);
    delete msg;
    compute_continue_();
    performance_stop_(perf_compute,__FILE__,__LINE__);
  }

  void p_compute_exit()
  {
    performance_start_(perf_compute,__FILE__,__LINE__);
    compute_exit_();
    performance_stop_(perf_compute,__FILE__,__LINE__);
  }
  void r_compute_exit(CkReductionMsg * msg)
  {
    performance_start_(perf_compute,__FILE__,__LINE__);
    delete msg;
    compute_exit_();
    performance_stop_(perf_compute,__FILE__,__LINE__);
    
  }

  /// Return the currently active Method
  int index_method() const throw()
  { return index_method_; }

  /// Return the currently-active Method
  Method * method () throw();

  /// Start a new solver
  void push_solver(int index_solver) throw()
  {  index_solver_.push_back(index_solver); }

  /// Return from a solver
  int pop_solver() throw()
  {
    int index = index_solver();
    ASSERT ("Block::pop_solver",
	    "Trying to pop element off of empty Block::index_solver_ stack",
	    index_solver_.size() > 0);
    index_solver_.resize(index_solver_.size()-1);
    return index;
  }

  /// Return the index of the current solver
  int index_solver() const throw()
  {
    ASSERT1("Block::index_solver()","%s index_solver_[] stack is empty",
	   name().c_str(),
	   (index_solver_.size() > 0));
    return index_solver_[index_solver_.size()-1];
  }
  
  /// Return the currently-active Solver
  Solver * solver () throw();

protected: // methods
  
  /// Enter control compute phase
  void compute_enter_();
  /// Initiate computing the sequence of Methods
  void compute_begin_();
  /// Initiate computing the next Method in the sequence
  void compute_next_();
  /// Return after performing any Refresh operations
  void compute_continue_();
  /// Cleanup after all Methods have been applied
  void compute_end_();
  /// Exit control compute phase
  void compute_exit_();

public: // methods

  /// Prepare to call compute_next_() after computing (used to
  /// synchronize between methods) Must be called at end of Method
  void compute_done();

  //--------------------------------------------------
  // OUTPUT
  //--------------------------------------------------

  void p_output_enter()
  {
    performance_start_(perf_output);
    output_enter_();
    performance_stop_(perf_output);
  }
  void r_output_enter(CkReductionMsg * msg)
  {
    performance_start_(perf_output);
    delete msg;
    output_enter_();
    performance_stop_(perf_output);

  }

  void p_output_end();

  void p_output_exit()
  {
    performance_start_(perf_output);
    output_exit_();
    performance_stop_(perf_output);
    
  }
  void r_output_exit(CkReductionMsg * msg)
  {
    performance_start_(perf_output);
    delete msg;
    output_exit_();
    performance_stop_(perf_output);
  }

  /// Contribute block data to ith output object in the simulation
  void p_output_write (int index_output);

  /// Contribute block data to the Initial input object
  void p_output_read (int index_input = 0)
  {
    performance_start_(perf_output);
    INCOMPLETE("Block::p_output_read");
    performance_stop_(perf_output);
  }

protected:
  void output_enter_();
  void output_begin_();
  void output_exit_();
public:

  //--------------------------------------------------
  // ADAPT
  //--------------------------------------------------

  void p_adapt_enter() 
  {
    performance_start_(perf_adapt_apply);
    adapt_enter_();
    performance_stop_(perf_adapt_apply);
    performance_start_(perf_adapt_apply_sync);
  }
  void r_adapt_enter(CkReductionMsg * msg) 
  {
    performance_start_(perf_adapt_apply);
    adapt_enter_(); delete msg;
    performance_stop_(perf_adapt_apply);
    performance_start_(perf_adapt_apply_sync);
  }

  void p_adapt_next ()
  {
    performance_start_(perf_adapt_update);
    adapt_next_();
    performance_stop_(perf_adapt_update);
    performance_start_(perf_adapt_update_sync);
  }
  void r_adapt_next (CkReductionMsg * msg)
  {
    performance_start_(perf_adapt_update);
    delete msg;    
    adapt_next_();
    performance_stop_(perf_adapt_update);
    performance_start_(perf_adapt_update_sync);
  }

  void p_adapt_called() 
  {
    performance_start_(perf_adapt_notify);
    adapt_called_();
    performance_stop_(perf_adapt_notify);
    performance_start_(perf_adapt_notify_sync);
  }
  void r_adapt_called(CkReductionMsg * msg) 
  {
    performance_start_(perf_adapt_notify);
    delete msg;    
    adapt_called_();
    performance_stop_(perf_adapt_notify);
    performance_start_(perf_adapt_notify_sync);
  }

  void p_adapt_end ()  
  {
    performance_start_(perf_adapt_end);
    adapt_end_();
    performance_stop_(perf_adapt_end);
    performance_start_(perf_adapt_end_sync);
  }
  void r_adapt_end (CkReductionMsg * msg)  
  {
    performance_start_(perf_adapt_end);
    delete msg;    
    adapt_end_();
    performance_stop_(perf_adapt_end);
    performance_start_(perf_adapt_end_sync);
  }

  void p_adapt_exit() 
  {
    performance_start_(perf_adapt_end);
    adapt_exit_();
    performance_stop_(perf_adapt_end);
    performance_start_(perf_adapt_end_sync);
  }
  void r_adapt_exit(CkReductionMsg * msg) 
  {
    performance_start_(perf_adapt_end);
    delete msg;
    adapt_exit_();
    performance_stop_(perf_adapt_end);
    performance_start_(perf_adapt_end_sync);
  }


  /// Parent tells child to delete itself
  void p_adapt_delete();
  void p_adapt_recv_level 
  (Index index_debug, 
   int ic3[3], 
   int if3[3],
   int level_now, int level_new);

  void p_adapt_recv_child (MsgCoarsen * msg);

  void adapt_recv (const int of3[3], const int ic3[3],
		   int level_face_new, int level_relative);

  void adapt_send_level();

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

  /// Syncronize before continuing with next callback
  void control_sync (int entry_point, int sync, int id = 0);

  /// synchronize with count other chares; count only needs to be supplied once
  void p_control_sync_count(int entry_point, int id, int count) 
  {
    performance_start_(perf_control);
    control_sync_count_(entry_point,id, count);
    performance_stop_(perf_control);
  }

protected:
  void control_sync_neighbor_(int entry_point, int id);
  void control_sync_face_(int entry_point, int id);
  void control_sync_count_(int entry_point, int id, int count);
public:

  //--------------------------------------------------
  // REFRESH
  //--------------------------------------------------

  void refresh_enter (int call, Refresh * refresh);

  /// Enter the refresh phase after synchronizing
  void p_refresh_continue ()
  {
    refresh_continue();
  }

  void refresh_continue();

  /// Exit the refresh phase after synchronizing
  void p_refresh_exit () 
  {
    performance_start_(perf_refresh_exit);
    refresh_exit_();
    performance_stop_(perf_refresh_exit);
    performance_start_(perf_refresh_exit_sync);
  }
  void r_refresh_exit (CkReductionMsg * msg) 
  {
    performance_start_(perf_refresh_exit);
    delete msg;    
    refresh_exit_();
    performance_stop_(perf_refresh_exit);
    performance_start_(perf_refresh_exit_sync);
  }
protected:
  void refresh_exit_ ();
public:

  /// Synchronize at start of refresh so that neighboring Blocks have
  /// all reached this point before exchanging data
  //  void p_refresh_sync()
  //  { refresh_sync_(0); }
protected:
  //  void refresh_sync_(int count);
  void refresh_load_faces_();

public:

  void p_refresh_store (MsgRefresh * msg);

  /// Get restricted data from child when it is deleted
  void p_refresh_child (int n, char a[],int ic3[3]);

protected:

  //--------------------------------------------------
  // REFRESH
  //--------------------------------------------------
  void refresh_begin_();

  /// Pack field face data into arrays and send to neighbors
  void refresh_load_field_faces_ (Refresh * refresh);

  /// Scatter particles in ghost zones to neighbors
  void refresh_load_particle_faces_ (Refresh * refresh);

  void refresh_load_field_face_
  (int refresh_type, Index index, int if3[3], int ic3[3]);
  void refresh_load_particle_face_
  (int refresh_type, Index index, int if3[3], int ic3[3]);

  //--------------------------------------------------
  // PARTICLES
  //--------------------------------------------------

  /// Create ParticleData objects for particle_array[] (with possibly
  /// duplicated pointers) and particle_list[] (with unique
  /// ParticleData objects) Also calls
  /// particle_determine_periodic_update_
  int particle_create_array_neighbors_
  (Refresh * refresh,
   ParticleData * particle_array[],
   ParticleData * particle_list[],
   Index index_list[]);

  /// Computes updates to positions (dpx[i],dpy[i],dpz[i]) for faces that
  /// cross periodic domain boundaries
  void particle_determine_periodic_update_ 
  (int * index_lower, int * index_upper, double * dpx, double * dpy, double * dpz);
  void particle_apply_periodic_update_
  ( int nl, ParticleData * particle_list[], Refresh * refresh);

  /// Scatter particles of given types in type_list, to appropriate
  /// particle_array ParticleData elements
  void particle_scatter_neighbors_
  (int npa, ParticleData * particle_array[],
   std::vector<int> & type_list, Particle particle_src);

  /// Scatter particles to appropriate partictle_list elements
  void particle_scatter_children_ (ParticleData * particle_list[],
				   Particle particle_src);

  /// Send particles in list to corresponding indices
  void particle_send_(int nl,Index index_list[], 
		      ParticleData * particle_list[]);

  /// Pack particle type data into arrays and send to neighbors
  void particle_load_faces_ (int npa, int * nl,
			     ParticleData * particle_list[],
			     ParticleData * particle_array[],
			     Index index_list[],
			     Refresh * refresh);

  //--------------------------------------------------
  // STOPPING
  //--------------------------------------------------

public:
  /// Entry method after begin_stopping() to call Simulation::r_stopping()
  void r_stopping_compute_timestep(CkReductionMsg * msg);

  /// Enter the stopping phase
  void p_stopping_enter () 
  {
    performance_start_(perf_stopping);
    stopping_enter_();
    performance_stop_(perf_stopping);
  }
  void r_stopping_enter (CkReductionMsg * msg) 
  {
    performance_start_(perf_stopping);
    delete msg;    
    stopping_enter_();
    performance_stop_(perf_stopping);
  }

  /// Quiescence before load balancing
  void p_stopping_balance();

  /// Manual call of LB
  void p_balance();

  /// Exit the stopping phase
  void p_stopping_exit () 
  {
    performance_start_(perf_stopping);
    stopping_exit_();
    performance_stop_(perf_stopping);
  }
  void r_stopping_exit (CkReductionMsg * msg) 
  {
    delete msg;    
    stopping_exit_();
    performance_stop_(perf_stopping);
  }
  
protected:

  void stopping_enter_();
  void stopping_begin_();
  void stopping_balance_();
  void stopping_exit_();

public:
  /// Exit the stopping phase to exit
  void p_exit () 
  {
    performance_start_(perf_exit);
    exit_();
    performance_stop_(perf_exit);
  }
  void r_exit (CkReductionMsg * msg) 
  {
    performance_start_(perf_exit);
    delete msg;    
    exit_();
    performance_stop_(perf_exit);
  }
protected:

  void exit_();

  //--------------------------------------------------
  // PERFORMANCE
  //--------------------------------------------------

protected:
  /// Start and stop measuring Block-based performance regions
  void performance_start_
  (int index_region, std::string file="", int line=0);
  void performance_stop_
  (int index_region, std::string file="", int line=0);

  //--------------------------------------------------
  // TESTING
  //--------------------------------------------------

  /// Check consistency between is_leaf_ and size of children_()
  void check_leaf_();

  /// Check if Block should have been deleted
  void check_delete_();

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

  /// Set Block's cycle
  virtual void set_cycle (int cycle) throw()
  { cycle_ = cycle;}

  /// Set Block's time
  virtual void set_time (double time) throw()
  { time_  = time; }

  /// Set Block's timestep
  virtual void set_dt (double dt) throw()
  { dt_  = dt; }

  /// Set Block's stopping criteria
  virtual void set_stop (double stop) throw()
  { stop_  = stop; }

  /// Initialize Block
  virtual void initialize () throw()
  {  }

  /// Return the local simulation object
  Simulation * simulation() const;

  /// Return the rank of the Simulation
  int rank() const;

  //  int count_neighbors() const;

  void ResumeFromSync();

  FieldFace * create_face
  (int if3[3], int ic3[3], bool lg3[3], int refresh_type,
   std::vector<int> & field_list);

protected: // functions

  /// Return the child adjacent to the given child in the direction of
  /// the given face
  void facing_child_(int jc3[3], const int ic3[3], const int if3[3]) const;

  /// Return whether the given face of the given child and its parent
  /// intersect
  bool child_is_on_face_(const int ic3[3],const int if3[3]) const;

  /// Return the face of the parent corresponding to the given face
  /// of the given child.
  bool parent_face_(int ipf3[3],const int if3[3], const int ic3[3]) const;

  void check_child_(const int ic3[3], const char * msg, 
		    const char * file, int line) const 
  {
    ASSERT5 (msg, "child %d %d %d out of range in file %s line %d",
	     ic3[0],ic3[1],ic3[2],file,line,
	     0 <= ic3[0] && ic3[0] <= 1 &&
	     0 <= ic3[1] && ic3[1] <= 1 &&
	     0 <= ic3[2] && ic3[2] <= 1);
  }

  void check_face_(const int if3[3], const char * msg,
		   const char * file, int line) const 
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

  /// Allocate and copy in attributes from give Block
  void copy_(const Block & block) throw();

  /// Return the (lower) indices of the Block in the level, 
  /// and the number of indices

  /// Update face_level_[] for given child in refined Block
  void refine_face_level_update_ (Index index_child);

  /// Update face_level_[] for coarsened Block
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

  /// Boundary is a boundary face
  bool is_boundary_face_(int of3[3],
			 bool boundary[3][2],
			 bool periodic[3][2],
			 bool update[3][2]) const 
  {

    bool skip = false;
    for (int axis=0; axis<3; axis++) {
      if (of3[axis] != 0) {
	int face=(of3[axis]+1)/2;
	if ( (! periodic[axis][face]) && 
	     update[axis][face] &&
	     boundary[axis][face]) skip = true;
      }
    }
    return skip;
  }

  /// Set the current refresh object
  void set_refresh (Refresh * refresh) 
  {  refresh_.push_back(*refresh); };

  /// Return the currently-active Refresh object
  Refresh * refresh () throw()
  {  return &refresh_.back();  }


  void print () const;
  
protected: // attributes

  /// Whether data exists
  /// Mesh Data that this Block controls
  Data * data_;

  /// Data for holding restricted child data 
  Data * child_data_;

  //--------------------------------------------------

  /// Index of this Block in the octree forest
  Index index_;

  /// Desired level for the next cycle
  int level_next_;

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

  /// Index of current initialization routine
  int index_initial_;

  /// MESH REFINEMENT

  /// list of child nodes
  std::vector<Index> children_;

  /// Synchronization counter for coarsening
  Sync sync_coarsen_;

  /// Synchronization counters for p_control_sync
  int  count_sync_[PHASE_COUNT];
  int  max_sync_[PHASE_COUNT];

  /// current level of neighbors along each face
  std::vector<int> face_level_curr_;

  /// new level of neighbors along each face
  std::vector<int> face_level_next_;

  /// current level of neighbors accumulated from children that can coarsen
  std::vector<int> child_face_level_curr_;

  /// new level of neighbors accumulated from children that can coarsen
  std::vector<int> child_face_level_next_;

  /// Can coarsen only if all children can coarsen
  int count_coarsen_;

  /// Number of adapt steps in the adapt phase
  int adapt_step_;

  /// Current adapt value for the block
  int adapt_;

  /// whether Block has been coarsened and should be deleted
  bool coarsened_;

  /// Whether Block is marked for deletion
  bool delete_;

  /// Whether Block is a leaf node during adapt phase (stored not
  /// computed to avoid race condition bug #30)
  bool is_leaf_;

  /// Age of the Block in cycles (for OutputImage)
  int age_;

  /// Last face level received from given face
  std::vector<int> face_level_last_;

  /// String for storing bit ID name
  mutable std::string name_;

  /// Index of currently-active Method
  int index_method_;

  /// Stack of currently active solvers
  std::vector<int> index_solver_;
  
  /// Refresh object associated with current refresh operation
  /// (Not a pointer since must be one per Block for synchronization counters)
  std::vector<Refresh> refresh_;

};

#endif /* COMM_BLOCK_HPP */

