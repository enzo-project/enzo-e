// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:01:51 PST 2011
/// @brief    [\ref Enzo] Declaration of the EnzoBlock class

#ifndef ENZO_ENZO_BLOCK_HPP
#define ENZO_ENZO_BLOCK_HPP

//----------------------------------------------------------------------

// #define TRACE_BLOCK

#include "charm_enzo.hpp"

class EnzoBlock : public CBase_EnzoBlock

{

  /// @class    EnzoBlock
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] An EnzoBlock is a Block with Enzo data

  friend class IoEnzoBlock;

  friend class EnzoSimulation;
  friend class EnzoTimestep;
  friend class EnzoTimestepPpm;
  friend class EnzoTimestepPpml;
  friend class EnzoInitialImplosion2;
  friend class EnzoInitialSedovArray2;

public: // interface

  /// Initialize the EnzoBlock chare array

  EnzoBlock ( process_type ip_source, MsgType msg_type );
  /// Initialize EnzoBlock using MsgRefine returned by creating process
  void set_msg_refine(MsgRefine * msg);
  void set_msg_check(EnzoMsgCheck * msg);

  /// Initialize an empty EnzoBlock
  EnzoBlock()
    :  CBase_EnzoBlock(),
       dt(0.0),
       redshift(0.0)
  {
    performance_start_(perf_block);

    for (int i=0; i<MAX_DIMENSION; i++) {
      GridLeftEdge[i] = 0;
      GridDimension[i] = 0;
      GridStartIndex[i] = 0;
      GridEndIndex[i] = 0;
      CellWidth[i] = 0.0;
    }
    performance_stop_(perf_block);
  }

  /// Charm++ Migration constructor
  EnzoBlock (CkMigrateMessage *m);

  /// Pack / unpack the EnzoBlock in a CHARM++ program
  void pup(PUP::er &p);

  /// Destructor
  virtual ~EnzoBlock();

  //--------------------------------------------------
  // Charm++ virtual
  //--------------------------------------------------

  virtual const CProxy_Block proxy_array() const
  { return thisProxy; }

  virtual const CProxyElement_Block proxy_element() const
  { return thisProxy[thisIndex]; }

  /// Write attributes, e.g. to stdout for debugging
  void write(FILE *fp=stdout) throw ();

  //--------------------------------------------------
  // Nested grid initialization.
  //--------------------------------------------------

  /// Check if this block should create child blocks during initialization.
  bool spawn_child_blocks() throw();

  /// Create child blocks.
  virtual void create_initial_child_blocks();
  void instantiate_children() throw();

  //----------------------------------------------------------------------
  // Original Enzo functions
  //----------------------------------------------------------------------

  /// Set EnzoBlock's dt (overloaded to update EnzoBlock::dt)
  virtual void set_dt (double dt) throw();

  /// Set EnzoBlock's time (overloaded to update current time)
  virtual void set_time (double time) throw();

  /// Set EnzoBlock's stopping criteria
  void set_stop (bool stop) throw();

  /// Initialize EnzoBlock
  virtual void initialize () throw();

public: /// entry methods

  /// Perform the necessary reductions
  CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs);

  /// Compute sum, min, and max of g values for EnzoMethodTurbulence
  void r_method_turbulence_end(CkReductionMsg *msg);

  void p_initial_hdf5_recv(MsgInitial * msg_initial);

  /// TEMP
  double timestep() { return dt; }

  //--------------------------------------------------

  // EnzoMethodBalance
  void p_method_balance_migrate();
  void p_method_balance_done();

  /// Synchronize after potential solve and before accelerations
  void p_method_gravity_continue();

  /// Synchronize for refresh
  void p_method_gravity_end();

  // EnzoMethodInference

  /// Merge inference array creation masks from children
  void p_method_infer_merge_masks (int n, char * mask, int ic3[3]);
  /// Count inference arrays
  void p_method_infer_count_arrays (int count);
  /// Request from inference array for field data
  void p_method_infer_request_data (int il3[3]);
  /// Update blocks with inference results
  void p_method_infer_update(int n, char * buffer, int il3[3]);
  /// Exit EnzoMethodInference
  void p_method_infer_exit();

  //--------------------------------------------------
  /// Checkpoint
  //--------------------------------------------------

  /// Call to Block array to self-identify as "first" when writing
  /// checkpoint files based on Ordering object
  void p_check_write_first(int num_files, std::string ordering, std::string);

  /// Call to single Block to return data for checkpoint
  void p_check_write_next(int num_files, std::string ordering);

  /// Exit EnzoMethodCheck
  void p_check_done();

  /// Initialize restart data in existing Block (level == 0)
  void p_restart_set_data(EnzoMsgCheck * );

  /// Refine to create the specified child in this block
  void p_restart_refine (int ic3[3], int io_reader, int ip);

  /// Exit restart
  void p_restart_done();

  //--------------------------------------------------

  /// Synchronize for accumulate refresh before adding the
  /// sink fields.
  void p_method_accretion_end();

  // ------------------------------------------------
  
  /// EnzoSolverCg entry method: DOT ==> refresh P
  void r_solver_cg_loop_0a (CkReductionMsg * msg);

  /// EnzoSolverCg entry method: ==> refresh P
  void r_solver_cg_loop_0b (CkReductionMsg * msg);

  /// EnzoSolverCg entry method: DOT(R,R) after shift
  void r_solver_cg_shift_1 (CkReductionMsg * msg);

  /// EnzoSolverCg entry method
  void p_solver_cg_loop_2 ();

  /// EnzoSolverCg entry method: DOT(P,AP)
  void r_solver_cg_loop_3 (CkReductionMsg * msg);

  /// EnzoSolverCg entry method: DOT(R,R)
  void r_solver_cg_loop_5 (CkReductionMsg * msg);

  /// EnzoSolverCg entry method:
  /// perform the necessary reductions for shift
  CkReductionMsg * r_solver_cg_shift(int n, CkReductionMsg ** msgs);

  void p_solver_cg_matvec();

  //--------------------------------------------------

  /// EnzoSolverBiCGStab entry method: SUM(B) and COUNT(B)
  void r_solver_bicgstab_start_1(CkReductionMsg* msg);

  /// EnzoSolverBiCGStab entry method: DOT(R,R)
  void r_solver_bicgstab_start_3(CkReductionMsg* msg);

  /// EnzoSolverBiCGStab entry method: return from preconditioner
  void p_solver_bicgstab_loop_2();

  /// EnzoSolverBiCGStab entry method: refresh Y
  void p_solver_bicgstab_loop_3();

  /// EnzoSolverBiCGStab entry method: DOT(V,R0), SUM(Y) and SUM(V)
  void r_solver_bicgstab_loop_5(CkReductionMsg* msg);

  /// EnzoSolverBiCGStab entry method: return from preconditioner
  void p_solver_bicgstab_loop_8();

  /// EnzoSolverBiCGStab entry method: refresh Y
  void p_solver_bicgstab_loop_9();

  /// EnzoSolverBiCGStab entry method: DOT(U,U), DOT(U,Q), SUM(Y) and SUM(U)
  void r_solver_bicgstab_loop_11(CkReductionMsg* msg);

  /// EnzoSolverBiCGStab entry method: DOT(R,R) and DOT(R,R0)
  void r_solver_bicgstab_loop_13(CkReductionMsg* msg);

  /// EnzoSolverBiCGStab entry method: ITER++
  void r_solver_bicgstab_loop_15(CkReductionMsg* msg);

  void p_dot_recv_parent  (int n, long double * dot_block,
			   std::vector<int> is_array,
			   int i_function, int iter);
  void p_dot_recv_children(int n, long double * dot_block,
			   std::vector<int> is_array,
			   int i_function);

/// EnzoSolverDd

  void p_solver_dd_restrict_recv(FieldMsg * msg);
  void p_solver_dd_prolong_recv(FieldMsg * msg);
  void solver_dd_prolong_recv(FieldMsg * msg);
  void p_solver_dd_solve_coarse();
  void p_solver_dd_solve_domain();
  void p_solver_dd_last_smooth();
  void r_solver_dd_barrier(CkReductionMsg* msg);
  void r_solver_dd_end(CkReductionMsg* msg);

  // EnzoSolverJacobi

  void p_solver_jacobi_continue();

  // EnzoSolverMg0

  void r_solver_mg0_begin_solve(CkReductionMsg* msg);
  void p_solver_mg0_restrict();
  void p_solver_mg0_solve_coarse();
  void p_solver_mg0_post_smooth();
  void p_solver_mg0_last_smooth();
  void r_solver_mg0_barrier(CkReductionMsg* msg);
  void p_solver_mg0_prolong_recv(FieldMsg * msg);
  void solver_mg0_prolong_recv(FieldMsg * msg);
  void p_solver_mg0_restrict_recv(FieldMsg * msg);

  // EnzoMethodFeedbackSTARSS
  void p_method_feedback_starss_end();

  //EnzoMethodM1Closure
  void p_method_m1_closure_solve_transport_eqn();
  void p_method_m1_closure_set_global_averages(CkReductionMsg * msg);

  virtual void print() const {
    FILE *fp = fopen ((std::string("EB-")+name_).c_str(),"a");
    fprintf (fp,"PRINT_ENZO_BLOCK name = %s\n",name().c_str());
    fprintf (fp,"PRINT_ENZO_BLOCK dt = %g\n",dt);
    fprintf (fp,"PRINT_ENZO_BLOCK redshift = %g\n",redshift);
    fprintf (fp,"PRINT_ENZO_BLOCK GridLeftEdge[] = %g %g %g\n",GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);
    fprintf (fp,"PRINT_ENZO_BLOCK GridDimension[] = %d %d %d\n",GridDimension[0],GridDimension[1],GridDimension[2]);
    fprintf (fp,"PRINT_ENZO_BLOCK GridStartIndex[] = %d %d %d\n",GridStartIndex[0],GridStartIndex[1],GridStartIndex[2]);
    fprintf (fp,"PRINT_ENZO_BLOCK GridEndIndex[] = %d %d %d\n",GridEndIndex[0],GridEndIndex[1],GridEndIndex[2]);
    fprintf (fp,"PRINT_ENZO_BLOCK CellWidth[] = %g %g %g\n",CellWidth[0],CellWidth[1],CellWidth[2]);
    Block::print(fp);
    fclose (fp);
  }

protected: // methods

  /// Create EnzoMsgCheck, returning file index
  int create_msg_check_
  ( EnzoMsgCheck ** msg_check, int num_files, std::string ordering,
    std::string name_dir = "", bool * is_first = nullptr);

  /// Initialize restart data in Block
  void restart_set_data_(EnzoMsgCheck * );

  /// Create a DataMsg object for this block
  DataMsg *create_data_msg_();

protected: // attributes

public: // attributes (YIKES!)

  union {
    enzo_float dt;
    enzo_float dtFixed;
  };

  /// Cosmological redshift for the current cycle
  enzo_float redshift;

  /// starting pos (active problem space)
  enzo_float GridLeftEdge[MAX_DIMENSION];
  /// total dimensions of all grids
  int GridDimension[MAX_DIMENSION];
  /// starting index of the active region
  int GridStartIndex[MAX_DIMENSION];
  /// stoping index of the active region
  int GridEndIndex[MAX_DIMENSION];
  enzo_float CellWidth[MAX_DIMENSION];

};

#endif /* ENZO_ENZO_BLOCK_HPP */

