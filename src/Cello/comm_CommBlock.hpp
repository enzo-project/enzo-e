// See LICENSE_CELLO file for license and copyright information

/// @file     comm_CommBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-28
/// @brief    [\ref Comm] Declaration of the CommBlock class

#ifndef COMM_COMM_BLOCK_HPP
#define COMM_COMM_BLOCK_HPP

class CommBlock {

  /// @class    CommBlock
  /// @ingroup  Comm
  /// @brief    [\ref Comm] 

public: // interface

  /// Constructor
  CommBlock() throw();

  /// Destructor
  ~CommBlock() throw();

  /// Copy constructor
  CommBlock(const CommBlock & CommBlock) throw();

  /// Assignment operator
  CommBlock & operator= (const CommBlock & CommBlock) throw();

  /// (X) Initialize block for the simulation.
  void p_initial();

  /// (X) Refresh ghost zones and apply boundary conditions
  void p_refresh() { refresh(); }

  /// (X) Apply the numerical methods on the block
  void p_compute(int cycle, double time, double dt);
  
  /// (X) Refresh a FieldFace
  void x_refresh(int n, char buffer[],int fx, int fy, int fz);

  /// (X) Contribute block data to ith output object in the simulation
  void p_write (int index_output);

  /// (X) Contribute block data to the Initial input object
  void p_read (int index_input = 0);

  /// (X) Entry function after prepare() to call Simulation::p_output()
  void p_output(CkReductionMsg * msg);


#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
  }
#endif
  
private: // functions


private: // attributes

  /// Block that this CommBlock is associated with
  Block * block_;
  // NOTE: change pup() function whenever attributes change

};

#endif /* COMM_COMM_BLOCK_HPP */

