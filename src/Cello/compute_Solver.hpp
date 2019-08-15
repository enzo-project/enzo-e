// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Solver.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Compute] Declaration for the Solver class

#ifndef COMPUTE_SOLVER_HPP
#define COMPUTE_SOLVER_HPP

#include <cstring>

class Refresh;
class Solver : public PUP::able 
{
  /// @class    Solver
  /// @ingroup  Compute
  /// @brief    [\ref Solver] Interface to a linear solver

public: // interface

  /// Create a new Solver
  Solver (std::string name,
	  std::string field_x,
	  std::string field_b,
	  int monitor_iter,
	  int restart_cycle,
	  int solve_type,
	  int min_level = 0,
	  int max_level = std::numeric_limits<int>::max()) throw();

  /// Create an uninitialized Solver
  Solver () throw();

  /// Charm++ PUP::able declarations
  PUPable_abstract(Solver);

  Solver (CkMigrateMessage *m)
    : PUP::able (m),
    name_(""),
    ix_(-1),ib_(-1),
    monitor_iter_(0),
    restart_cycle_(1),
    callback_(0),
    index_(0),
    min_level_(- std::numeric_limits<int>::max()),
    max_level_(  std::numeric_limits<int>::max()),
    id_sync_(0),
    solve_type_(solve_leaf)
#ifdef NEW_REFRESH
    ,ir_post_(-1)
#else /* ! NEW_REFRESH */
    ,refresh_list_()
#endif      
  { }

  /// Destructor
  virtual ~Solver() throw();
  
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    
    PUP::able::pup(p);

    p | name_;
    p | ix_;
    p | ib_;
    p | monitor_iter_;
    p | restart_cycle_;
    p | callback_;
    p | index_;
    p | min_level_;
    p | max_level_;
    p | id_sync_;
    p | solve_type_;
#ifdef NEW_REFRESH
    p | ir_post_;
#else    
    p | refresh_list_;
#endif    
  }

  Refresh * refresh(size_t index=0) ;

  void set_callback (int callback)
  { callback_ = callback; }

  void set_index (int index)
  { index_ = index; }

  int index() const
  { return index_; }

  void set_field_x (int ix)
  { ix_ = ix;  }
  
  void set_field_b (int ib)
  { ib_ = ib;  }

  void set_min_level (int min_level)
  { min_level_ = min_level; }

  int min_level()
  { return min_level_; }
  
  int max_level()
  { return max_level_; }

  void set_max_level (int max_level)
  { max_level_ = max_level; }

  void set_sync_id (int sync_id)
  { id_sync_ = sync_id; }
  
  /// Type of neighbor: level if min_level == max_level, else leaf
  int neighbor_type_() const throw() {
    int retval;
    switch (solve_type_) {
    case solve_level:
      retval = neighbor_level;
      break;
    case solve_leaf:
      retval = neighbor_leaf;
      break;
    case solve_tree:
      retval = neighbor_tree;
      break;
    default:
      retval = neighbor_unknown;
      break;
    }
    return retval;
  }

  /// Type of synchronization: sync_face if min_level == max_level,
  /// else sync_neighbor
  int sync_type_() const throw() {
    int retval;
    switch (solve_type_) {
    case solve_leaf:
      retval = sync_neighbor;
      break;
    case solve_level:
      retval = sync_face;
      break;
    case solve_tree:
      retval = sync_neighbor;
      break;
    default:
      retval = sync_unknown;
      break;
    }
    return retval;
  }

  /// Which subset of Blocks the solver is defined on; see enum solve_type
  /// for supported types
  int solve_type() const
  { return solve_type_; }

  /// Return the name of this solver
  std::string name () const
  { return name_; }

public: // virtual functions

  /// Solve the linear system Ax = b
  virtual void apply ( std::shared_ptr<Matrix> A, Block * block) throw() = 0;

  /// Return the type of this solver
  virtual std::string type () const = 0;

  /// Whether Block is active
  virtual bool is_active_(Block * block) const;

  /// Whether solution is defined on this Block
  virtual bool is_finest_(Block * block) const;
  
protected: // functions

  /// Initialize a solve
  void begin_(Block * block);
  
  /// Clean up after a solver is done and returning to its callback_
  void end_(Block * block);

  void monitor_output_(Block * block, int iter,
		       double rr0=0.0,
		       double rr_min=0.0, double rr=0.0, double rr_max=0.0,
		       bool final = false) throw();
  
#ifdef NEW_REFRESH

  /// Add a new refresh object
  int add_new_refresh_ ();

  /// Return the specified Refresh object
  Refresh & new_refresh(int ir);
  
  /// Return the index of the main post-refresh object
  int refresh_post_id() const;

  /// Return the main post-refresh object for the solver
  Refresh & refresh_post();

#else
  
  int add_refresh (int ghost_depth, 
		   int min_face_rank, 
		   int neighbor_type, 
		   int sync_type,
		   int sync_id);
#endif
  
  /// Perform vector copy X <- Y
  template <class T>
  void copy_ (T * X, const T * Y,
	      int mx, int my, int mz,
	      bool active = true) const throw()
  {
    if (! active ) return;
    const int m = mx*my*mz;
    std::memcpy(X,Y,m*sizeof(T));
  }

  int sync_id_() const throw()
  { return this->id_sync_; }

  bool reuse_solution_ (int cycle) const throw();

protected: // attributes

  /// Name of the solver
  std::string name_;

  /// Field id for solution
  int ix_;
  
  /// Field id for right-hand side
  int ib_;
  
  /// How often to write output
  int monitor_iter_;

  /// Whether to reuse the previous solution as the initial guess
  int restart_cycle_;

  /// Callback id
  int callback_;

  /// Index of this solver
  int index_;

  /// Minimum mesh level
  int min_level_;
  
  /// Maximum mesh level
  int max_level_;

  /// Sync id
  int id_sync_;

  /// Type of solver; see enum solve_type for supported types
  int solve_type_;

#ifdef NEW_REFRESH
  /// New Refresh id for after the solver
  int ir_post_;
#else
  ///  Refresh object
  std::vector<Refresh *> refresh_list_;

#endif  
};

#endif /* COMPUTE_SOLVER_HPP */
