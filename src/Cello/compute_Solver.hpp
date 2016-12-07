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
  Solver () throw()
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

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    
    PUP::able::pup(p);
    
    p | refresh_list_;
  }

  int add_refresh (int ghost_depth, 
		   int min_face_rank, 
		   int neighbor_type, 
		   int sync_type);

  Refresh * refresh(size_t index=0) ;

public: // virtual functions

  /// Solve the linear system Ax = b

  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw() = 0; 

protected: // attributes

  ///  Refresh object
  std::vector<Refresh *> refresh_list_;

};

#endif /* COMPUTE_SOLVER_HPP */
