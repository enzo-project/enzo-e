// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Method.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Problem] Declaration for the Method class

#ifndef PROBLEM_METHOD_HPP
#define PROBLEM_METHOD_HPP

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
  {}

  /// Charm++ PUP::able declarations
  PUPable_abstract(Method);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);
    p | refresh_;
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

  Refresh * refresh() { return refresh_; }

protected: // functions

  /// Perform vector copy X <- Y
  template <class T>
  void copy_ (T * X, const T * Y,
	      int mx, int my, int mz,
	      bool is_active) const throw()
  {
    if (! is_active ) return;
    const int m = mx*my*mz;
    for (int i=0; i<m; i++) X[i] = Y[i];
  }

protected: // attributes

  ///  Refresh object
  Refresh * refresh_;

};

#endif /* PROBLEM_METHOD_HPP */
