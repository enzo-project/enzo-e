// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Solver.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Compute] Declaration for the Solver class

#ifndef COMPUTE_SOLVER_HPP
#define COMPUTE_SOLVER_HPP

class Refresh;
class Solver : public PUP::able 
{
  /// @class    Solver
  /// @ingroup  Compute
  /// @brief    [\ref Solver] Interface to a linear solver

public: // interface

  /// Create a new Solver
  Solver (int monitor_iter) throw()
    : PUP::able(),
      monitor_iter_(monitor_iter)
  {}

  /// Create an uninitialized Solver
  Solver () throw()
    : PUP::able(),
      monitor_iter_(0)
  {}

  /// Destructor
  virtual ~Solver() throw()
  {
    for (size_t i=0; i<refresh_list_.size(); i++) {
      delete refresh_list_[i];
      refresh_list_[i] = 0;
    }
  }

  /// Charm++ PUP::able declarations
  PUPable_abstract(Solver);

  Solver (CkMigrateMessage *m)
    : PUP::able (m)
  { }
  
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    
    PUP::able::pup(p);
    
    p | refresh_list_;
    p | monitor_iter_;
    p | callback_;
    p | index_;
  }

  int add_refresh (int ghost_depth, 
		   int min_face_rank, 
		   int neighbor_type, 
		   int sync_type,
		   int id_in=0);

  Refresh * refresh(size_t index=0) ;

  void set_callback (int callback)
  { callback_ = callback; }

  void set_index (int index)
  { index_ = index; }
  
public: // virtual functions

  /// Solve the linear system Ax = b
  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw() = 0;

  /// Return the name of this solver
  virtual std::string name () const
  { return "UNKNOWN"; }

protected: // functions

  /// Initialize a solve
  void begin_(Block * block);
  
  /// Clean up after a solver is done and returning to its callback_
  void end_(Block * block);

  void monitor_output_(Block * block, int iter, double rr0,
		       double rr_min, double rr, double rr_max,
		       bool final = false) throw();
  
  /// Perform vector copy X <- Y
  template <class T>
  void copy_ (T * X, const T * Y,
	      int mx, int my, int mz,
	      bool active = true) const throw()
  {
    if (! active ) return;
    const int m = mx*my*mz;
    for (int i=0; i<m; i++) X[i] = Y[i];
  }

protected: // attributes

  ///  Refresh object
  std::vector<Refresh *> refresh_list_;

  /// How often to write output
  int monitor_iter_;

  /// Callback id
  int callback_;

  /// Index of this solver
  int index_;
};

#endif /* COMPUTE_SOLVER_HPP */
