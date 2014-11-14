// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 14:09:42 PDT 2010
/// @brief    [\ref Comm] Declaration of the CommBlock class
/// @todo     Move Enzo Gravity to EnzoBlock

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
// class Control;
class Factory;
class FieldDescr;
class FieldFace;
class ItFace;
class GroupProcess;
class Hierarchy;
class Method;
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
  CommBlock()  { printf ("CommBlock()\n"); };

  /// Initialize a migrated CommBlock
  CommBlock (CkMigrateMessage *m);

  /// CHARM pupper
  void pup(PUP::er &p);

  //----------------------------------------------------------------------
  // ACCESS METHODS
  //----------------------------------------------------------------------
  
  /// Return the Block associated with this CommBlock
  inline Block * block() throw()
  { return block_; };

  const FieldDescr * field_descr() throw();

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

  /// Return whether this CommBlock is a leaf in the octree forest
  bool is_leaf() const 
  { return is_leaf_; }

  /// Index of the CommBlock
  const Index & index() const 
  { return index_; }

  const int & face_level (const int if3[3]) const
  { return face_level_curr_[IF3(if3)]; }

  const int & face_level_next (const int if3[3]) const
  { return face_level_next_[IF3(if3)]; }

  const int & child_face_level (const int ic3[3], const int if3[3]) const
  { return child_face_level_curr_[ICF3(ic3,if3)]; }

  const int & child_face_level_next (const int ic3[3], const int if3[3]) const
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

  ItFace it_face(const int * ic3=0,
		 const int * if3=0) throw();

  //--------------------------------------------------
  // Charm++ virtual
  //--------------------------------------------------

  virtual const CProxy_CommBlock proxy_array() const 
  { return thisProxy; }

  virtual const CProxyElement_CommBlock proxy_element() const 
  { return thisProxy[thisIndex]; }

  //--------------------------------------------------
  // INITIAL
  //--------------------------------------------------

  void initial_exit_();
  void p_initial_exit()
  {      initial_exit_(); }
  void r_initial_exit(CkReductionMsg * msg)
  {      initial_exit_(); delete msg; }

  //--------------------------------------------------
  // COMPUTE
  //--------------------------------------------------

  void p_compute_enter()
  {      compute_enter_(); }
  void r_compute_enter(CkReductionMsg * msg)
  {      compute_enter_(); delete msg; }

  void p_compute_continue()
  {      compute_continue_(); }
  void r_compute_continue(CkReductionMsg * msg)
  {      compute_continue_(); delete msg; }

  void p_compute_exit()
  {      compute_exit_(); }
  void r_compute_exit(CkReductionMsg * msg)
  {      compute_exit_(); delete msg; }

  /// Set the currently-active Method.  Used to resume a Method's
  /// computation after a reduction
  void set_method (Method * method) throw()
  { method_ = method; }

  /// Return the currently-active Method
  Method * method () throw()
  { return method_; }

protected:
  void compute_enter_();
  void compute_continue_();
  void compute_begin_();
  void compute_exit_();
public:

  //--------------------------------------------------
  // ENZO GRAVITY CG: MOVE TO ENZO
  //--------------------------------------------------

  void p_enzo_matvec()
  {      enzo_matvec_(); }
  void r_enzo_matvec(CkReductionMsg * msg)
  {      enzo_matvec_(); delete msg; }

protected:
  void enzo_matvec_() 
  {
    printf ("%s:%d enzo_matvec_\n",__FILE__,__LINE__);
  }

public:


  //--------------------------------------------------
  // OUTPUT
  //--------------------------------------------------

  void p_output_enter()
  {      output_enter_(); }
  void r_output_enter(CkReductionMsg * msg)
  {      output_enter_(); delete msg;}

  void p_output_end();

  void p_output_exit()
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
  void r_adapt_enter(CkReductionMsg * msg) 
  {      adapt_enter_(); delete msg;}

  void p_adapt_next ()
  {      adapt_next_(); }
  void r_adapt_next (CkReductionMsg * msg)
  {      adapt_next_(); delete msg; }

  void p_adapt_called() 
  {      adapt_called_(); }
  void r_adapt_called(CkReductionMsg * msg) 
  {      adapt_called_(); delete msg; }

  void p_adapt_end ()  
  {      adapt_end_(); }
  void r_adapt_end (CkReductionMsg * msg)  
  {      adapt_end_(); delete msg;}

  void p_adapt_exit() 
  {      adapt_exit_(); }
  void r_adapt_exit(CkReductionMsg * msg) 
  {      adapt_exit_(); delete msg;}


  /// Parent tells child to delete itself
  void p_adapt_delete();
  void p_adapt_recv_level 
  (Index index_debug, int ic3[3], int if3[3], int level_now, int level_new);
  void p_adapt_recv_child
  (int ic3[3],int na, char * array, int nf, int * child_face_level);

  void adapt_recv_same   (const int of3[3],int level_face_new);
  void adapt_recv_coarse (const int of3[3],const int ic3[3],int level_face_new);
  void adapt_recv_fine   (const int of3[3],const int ic3[3],int level_face_new);
  void adapt_recv_recurse(const int if3[3],const int ic3[3],
			  int level_face,
			  int level_face_new,
			  Index index_send);

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

  /// Syncronize before continuing with next phase
  void control_next(int phase = phase_unknown, std::string sync = "");

  /// Syncronize before continuing with next phase
  void control_sync(int phase, std::string sync, bool next_phase, const char * name, int line);

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

  // Refresh ghost zones and apply boundary conditions
  void p_control_refresh_enter()  
  {      control_refresh_enter_(); }
  void r_control_refresh_enter(CkReductionMsg * msg)  
  {      control_refresh_enter_(); delete msg; }
  void control_refresh_enter_();

  // Exit the refresh phase after QD
  void p_control_refresh_exit () 
  {      control_refresh_exit_(); }
  void r_control_refresh_exit (CkReductionMsg * msg) 
  {      control_refresh_exit_(); delete msg;  }
  void control_refresh_exit_ ();

  /// Refresh a FieldFace in same, next-coarser, or next-finer level
  void x_refresh_send_face
  (int n, char buffer[],  int type_refresh, int if3[3], int ic3[3]) 
  {      refresh_store_face_(n,buffer,type_refresh,if3,ic3); }

  /// Get restricted data from child when it is deleted
  void x_refresh_child (int n, char buffer[],int ic3[3]);

  /// Whether the given face should be refreshed 
  /// (false if boundary and not periodic)
  bool do_refresh(int if3[3]) const;

protected:
  void refresh_begin_();
  void refresh_load_face_
  (int type_refresh, Index index, int if3[3], int ic3[3]);
  void refresh_store_face_
  (int n, char buffer[],  int type_refresh, int if3[3], int ic3[3]);
public:

  //--------------------------------------------------
  // STOPPING
  //--------------------------------------------------

  /// Entry method after begin_stopping() to call Simulation::r_stopping()
  void r_stopping_compute_timestep(CkReductionMsg * msg);

  /// Enter the stopping phase
  void p_stopping_enter () 
  {      stopping_enter_(); }
  void r_stopping_enter (CkReductionMsg * msg) 
  {      stopping_enter_(); delete msg;  }

  /// Enter the stopping phase
  void p_stopping_exit () 
  {      stopping_exit_(); }
  void r_stopping_exit (CkReductionMsg * msg) 
  {      stopping_exit_(); delete msg;  }

protected:
  void stopping_enter_();
  void stopping_begin_();
  void stopping_balance_();
  void stopping_exit_();
public:

  /// Exit the stopping phase to exit
  void p_exit () 
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

  /// Return the rank of the Simulation
  int rank() const;

  //  int count_neighbors() const;

  void ResumeFromSync();

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

  /// Current phase of the cycle
  int index_cycle_phase_;

  /// whether CommBlock has been coarsened and should be deleted
  bool coarsened_;

  /// Whether CommBlock is marked for deletion
  bool delete_;

  /// Whether CommBlock is a leaf node during adapt phase (stored not
  /// computed to avoid race condition bug #30)
  bool is_leaf_;

  /// Age of the CommBlock in cycles (for OutputImage)
  int age_;

  /// Last face level received from given face
  std::vector<int> face_level_last_;

  /// String for storing bit ID name
  std::string name_;

  /// Currently-active Method
  Method * method_;

  /// Phase after current refresh
  int refresh_phase_;

  /// Synchronization after current refresh
  std::string refresh_sync_;

};

#endif /* COMM_COMMBLOCK_HPP */

