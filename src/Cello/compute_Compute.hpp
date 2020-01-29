// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Compute.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Compute] Declaration for the Compute class

#ifndef COMPUTE_COMPUTE_HPP
#define COMPUTE_COMPUTE_HPP

class Compute : public PUP::able 
{
  /// @class    Compute
  /// @ingroup  Compute
  /// @brief    [\ref Compute] Interface to an application compute / analysis / visualization function.

public: // interface

  /// Create a new Compute
  Compute () throw()
  {
    set_history(0);
  }

  /// Destructor
  virtual ~Compute() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_abstract(Compute);

  Compute(CkMigrateMessage *m)
    : PUP::able(m),
    i_hist_(0)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);

    p | i_hist_;
  }

public: // virtual functions

  /// Perform the computation on the Block

  virtual void compute ( Block * block) throw() = 0; 

  /// Return the name of this Compute
  ///  used for associating derived fields when relevant
  virtual std::string name () throw () {
    return std::string();
  };

  /// Return / set field history to use in computation

  virtual int  get_history(int i_hist) {return i_hist_;};
  virtual void set_history(int i_hist) {i_hist_ = i_hist;};

protected:

  int i_hist_;

};

#endif /* COMPUTE_COMPUTE_HPP */
