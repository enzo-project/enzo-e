// See LICENSE_CELLO file for license and copyright information

/// @file     problem_InitialTrace.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2015-11-06 23:00:30
/// @brief    [\ref Problem] Declaration for the InitialTrace component

#ifndef PROBLEM_INITIAL_TRACE_HPP
#define PROBLEM_INITIAL_TRACE_HPP

class Hierarchy;

class InitialTrace : Initial
{

  /// @class    InitialTrace
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Initialize trace particles

public: // interface

  /// empty constructor for charm++ pup()
  InitialTrace() throw() {}

  /// Destructor
  virtual ~InitialTrace() throw()
  {} ;

  /// CHARM++ PUP::able declaration
  PUPable_decl(InitialTrace);

  /// CHARM++ migration constructor for PUP::able
  InitialTrace (CkMigrateMessage *m) : Initial(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);


public: // virtual functions

  /// InitialTraceize a Block
  virtual void enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const Hierarchy  * hierarchy
    ) throw();

protected: // functions


protected: // attributes

};

#endif /* PROBLEM_INITIAL_TRACE_HPP */
