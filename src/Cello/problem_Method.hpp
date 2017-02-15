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
  Method (double courant = 1.0) throw() : schedule_(NULL), courant_(courant) {}

  /// Destructor
  virtual ~Method() throw();

  /// Charm++ PUP::able declarations
  PUPable_abstract(Method);
  
  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual functions

  /// Apply the method to advance a block one timestep 

  virtual void compute ( Block * block) throw() = 0; 

  /// Return the name of this Method
  virtual std::string name () throw () = 0;

  /// Compute maximum timestep for this method
  virtual double timestep (Block * block) const throw() 
  { return std::numeric_limits<double>::max(); }

  /// Resume computation after a reduction
  virtual void compute_resume ( Block * block,
				CkReductionMsg * msg) throw()
  {
    /* This function intentionally empty */
  }

  int add_refresh (int ghost_depth, 
		   int min_face_rank, 
		   int neighbor_type, 
		   int sync_type,
		   int id=0)
  {
    int index=refresh_list_.size();
    refresh_list_.resize(index+1);
    refresh_list_[index] = new Refresh 
      (ghost_depth,min_face_rank,neighbor_type,sync_type,id,true);
    return index;
  }

  Refresh * refresh(size_t index=0) 
  {
    return (index < refresh_list_.size()) ? refresh_list_[index] : NULL;
  }

  /// Return the Schedule object pointer
  Schedule * schedule() throw() 
  { return schedule_; };

  /// Set schedule
  void set_schedule (Schedule * schedule) throw();

  double courant() const throw ()
  { return courant_; }

  void set_courant(double courant) throw ()
  { courant_ = courant; }

protected: // functions

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

  ///  Refresh object
  std::vector<Refresh *> refresh_list_;


  /// Schedule object, if any (default is every cycle)
  Schedule * schedule_;

  /// Courant condition for the Method
  double courant_;

};

#endif /* PROBLEM_METHOD_HPP */
