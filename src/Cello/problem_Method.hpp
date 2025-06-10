// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 13 11:11:47 PDT 2009
/// @brief    [\ref Problem] Declaration for the Method class

#ifndef PROBLEM_METHOD_HPP
#define PROBLEM_METHOD_HPP

class Refresh;
class Schedule;

class Method : public PUP::able
{
  /// @class    Method
  /// @ingroup  Method
  /// @brief    [\ref Method] Interface to an application method / analysis / visualization function.

public: // interface

  /// Create a new Method
  Method (double courant = 1.0) throw();

  /// Charm++ PUP::able declarations
  PUPable_abstract(Method);

  Method (CkMigrateMessage *m)
    : PUP::able(m),
      schedule_(NULL),
      courant_(1.0),
      ir_post_(-1),
      neighbor_type_(neighbor_leaf),
      max_supercycle_(1),
      index_method_(-1)
  { }

  /// Destructor
  virtual ~Method() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual functions

  /// Apply the method to advance a block one timestep
  ///
  /// The `Block::compute_done()` method MUST be invoked on all `block`s passed
  /// to this member function by the end of the control flow that this function
  /// launches.
  /// - In simple cases, that should be done just before this function returns.
  /// - The placement will varies in more complex cases. For example, if this
  ///   function invokes a single reduction, the call to
  ///   `Block::compute_done()` should be performed after completing the
  ///   reduction (e.g. in the `compute_resume` member function)
  virtual void compute ( Block * block) throw() = 0;

  /// Return the name of this Method
  virtual std::string name () throw () = 0;

  /// Compute maximum timestep for this method
  ///
  /// The default implementation returns the maximum finite value of `double`
  virtual double timestep (Block * block) throw()
  { return std::numeric_limits<double>::max(); }

  /// Resume computation after a reduction
  ///
  /// This member function only typically needs to be implemented by Method
  /// classes that employ reductions.
  virtual void compute_resume ( Block * block,
				CkReductionMsg * msg) throw()
  {
    /* This function intentionally empty */
  }

  /// Add a new refresh object
  int add_refresh_ (int neighbor_type = neighbor_leaf);

  /// Return the index for the main post-refresh object
  int refresh_id_post() const;

  /// Return the Schedule object pointer
  Schedule * schedule() throw()
  { return schedule_; };

  /// Set schedule
  void set_schedule (Schedule * schedule) throw();

  /// Query the associated courant factor
  double courant() const throw ()
  { return courant_; }

  void set_courant(double courant) throw ()
  { courant_ = courant; }

  /// Set maximum cycles for super-cycling
  void set_max_supercycle (int max_supercycle)
  { max_supercycle_ = max_supercycle; }

  /// Return maximum cycles for super-cycling
  int max_supercycle () const
  { return max_supercycle_; }

  /// Whether super-cycling is enabled
  bool is_supercycle() const
  { return (max_supercycle_ > 1); }

  /// Define a field and its two previously saved values for use
  /// in super-cycling, and return id_super identifying the values
  int super_define_fields_
  (std::string field,
   std::string field_curr,
   std::string field_prev);

  /// Save (curr = latest) or shift (prev = curr) the computed
  /// field(s) to use for extrapolating to get approximation in
  /// non-active cycles when super-cycling
  void super_save_fields_(Block * );
  void super_shift_fields_(Block * );

  /// Extrapolate using saved fields to get approximation at
  /// given time in non-active cycles when super-cycling
  void super_extrapolate_fields_(Block * block, double time );

  void super_update_time_(Block * block, double time);

void set_index(int index)
  { index_method_ = index; }

  int index() const
  { return index_method_; }

protected: // functions

  /// Whether this is a "solve-cycle" when supercycling
  bool is_solve_cycle_(Block * block);

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

public: // attributes (static)

  static double courant_global;

protected: // attributes

  /// Schedule object, if any (default is every cycle)
  Schedule * schedule_;

  /// Courant condition for the Method
  double courant_;

  /// Index for main refresh after Method is called
  int ir_post_;

  /// Default refresh type
  int neighbor_type_;

  /// Maximum allowed super-cycling
  int max_supercycle_;

  /// Index of this Method in Problem
  int index_method_;

  /// Field ID's used in super-cycling
  std::vector<int> super_field_;
  std::vector<int> super_field_curr_;
  std::vector<int> super_field_prev_;
  /// Scalars for times of saved fields
  int is_time_curr_;
  int is_time_prev_;

};

#endif /* PROBLEM_METHOD_HPP */
