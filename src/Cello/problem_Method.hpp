// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Problem] Declaration for the Method class

#ifndef PROBLEM_METHOD_HPP
#define PROBLEM_METHOD_HPP

class Refresh;

class Method : public PUP::able 
{
  /// @class    Method
  /// @ingroup  Method
  /// @brief    [\ref Method] Interface to an application method / analysis / visualization function.

public: // interface

  /// Create a new Method
  Method () throw()
  {}

  /// Destructor
  virtual ~Method() throw()
  {
    for (size_t i=0; i<refresh_list_.size(); i++) {
      delete refresh_list_[i];
      refresh_list_[i] = 0;
    }
  }

  /// Charm++ PUP::able declarations
  PUPable_abstract(Method);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);
    p | refresh_list_;
  }

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

  void set_num_refresh(size_t n)
  {
    // Reserve if needed
    if (n > refresh_list_.size()) {
      // if reserving, initialize only new elements
      int i1 = refresh_list_.size();
      int i2 = n;
      refresh_list_.resize(n);
      for (int i=i1; i < i2; i++) {
	refresh_list_[i] = new Refresh;
      }
    }
  }

  Refresh * refresh(size_t index=0) 
  {
    return (index < refresh_list_.size()) ? refresh_list_[index] : NULL;
  }

protected: // functions

  /// Perform vector copy X <- Y
  template <class T>
  void copy_ (T * X, const T * Y,
	      int mx, int my, int mz,
	      bool active) const throw()
  {
    if (! active ) return;
    const int m = mx*my*mz;
    for (int i=0; i<m; i++) X[i] = Y[i];
  }

protected: // attributes

  ///  Refresh object
  std::vector<Refresh *> refresh_list_;

};

#endif /* PROBLEM_METHOD_HPP */
