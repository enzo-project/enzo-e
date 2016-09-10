// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Initial.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Problem] Declaration for the Initial component

#ifndef METHOD_INITIAL_HPP
#define METHOD_INITIAL_HPP

class Hierarchy;

class Initial : public PUP::able 
{

  /// @class    Initial
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate an initial conditions generator

public: // interface

  /// empty constructor for charm++ pup()
  Initial() throw()
  : cycle_(0), time_(0.0) {}

  /// Create a new Initial
  Initial(int cycle, double time) throw()
    : cycle_(cycle), time_(time)
  {};

  /// Destructor
  virtual ~Initial() throw()
  {} ;

  /// CHARM++ PUP::able declaration
  PUPable_decl(Initial);

  /// CHARM++ migration constructor for PUP::able
  Initial (CkMigrateMessage *m)
    : PUP::able(m),
      cycle_(0),
      time_(0.0)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Initial time
  double time() const throw() { return time_; }

  /// Initial cycle
  int cycle() const throw() { return cycle_; }

public: // virtual functions

  /// Initialize a Block
  virtual void enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw()
  { enforce_block_(block,field_descr,particle_descr,hierarchy); }

  /// Return whether enforce() expects block != NULL
  virtual bool expects_blocks_allocated() const throw()
  { return true; }

protected: // functions

  /// Loop over enforcing initial conditions Field data in the Block
  void enforce_block_
  ( Block            * block, 
    const FieldDescr * field_descr,  
    const ParticleDescr * particle_descr,  
    const Hierarchy  * hierarchy
    ) throw();

protected: // attributes

  /// Initial cycle number
  int cycle_;

  /// Initial time
  double time_;

};

#endif /* METHOD_INITIAL_HPP */
