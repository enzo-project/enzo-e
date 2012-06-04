// See LICENSE_CELLO file for license and copyright information

/// @file     method_Initial.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Method] Declaration for the Initial component

#ifndef METHOD_INITIAL_HPP
#define METHOD_INITIAL_HPP

class Initial {

  /// @class    Initial
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate an initial conditions generator

public: // interface

  /// Create a new Initial
  Initial(int cycle, double time) throw()
    : cycle_(cycle), time_(time)
  {};

  /// Destructor
  virtual ~Initial() throw()
  {} ;

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    p | cycle_;
    p | time_;

  }
#endif


  /// Initial time
  double time() const throw() { return time_; }

  /// Initial cycle
  int cycle() const throw() { return cycle_; }

public: // virtual functions

#ifdef NEW_INITIAL

  /// Initialize an entire simulation
  virtual void enforce ( Simulation * simulation ) throw()
  { enforce_(simulation); }

  /// Initialize a mesh Hierarchy
  virtual void enforce
  ( Hierarchy * hierarchy, 
    const FieldDescr * field_descr  ) throw()
  { enforce_(hierarchy,field_descr); }

  /// Initialize a Patch
  virtual void enforce
  ( Patch * patch, 
    const FieldDescr * field_descr,
    const  Hierarchy * hierarchy
    ) throw()
  { enforce_(patch,field_descr,hierarchy); }

  /// Initialize a Block
  virtual void enforce
  ( Block * block, 
    const FieldDescr * field_descr,
    const Hierarchy * hierarchy
    ) throw()
  { enforce_(block,field_descr,hierarchy); }

#else /* not NEW_INITIAL */

  /// Enforce initial conditions on the given Block or Hierarchy
  virtual void enforce (Block * block,
			const FieldDescr * field_descr,
			const Hierarchy * hierarchy) throw() = 0;
#endif

  /// Return whether enforce() expects block != NULL
  virtual bool expects_blocks_allocated() const throw()
  { return true; }

protected: // functions

#ifdef NEW_INITIAL

  /// "Loop" over enforcing initial conditions on the Hierarchy
  void enforce_ (Simulation * simulation ) throw();

  /// Loop over enforcing initial conditions Patches in the Hierarchy
  void enforce_
  ( Hierarchy * hierarchy, 
    const FieldDescr * field_descr  ) throw();

  /// Loop over enforcing initial conditions Blocks in the Patch
  void enforce_
  ( Patch * patch, 
    const FieldDescr * field_descr,
    const Hierarchy * hierarchy
    ) throw();

  /// Loop over enforcing initial conditions Field data in the Block
  void enforce_
  ( Block * block, 
    const FieldDescr * field_descr,  
    const Hierarchy * hierarchy
    ) throw();


#endif

protected: // attributes

  /// Initial cycle number
  int cycle_;

  /// Initial time
  double time_;

};

#endif /* METHOD_INITIAL_HPP */
