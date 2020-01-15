// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodCloseFiles.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-01-14
/// @brief    [\ref Problem] Declaration of the MethodCloseFiles class

#ifndef PROBLEM_METHOD_CLOSE_FILES_HPP
#define PROBLEM_METHOD_CLOSE_FILES_HPP

class MethodCloseFiles : public Method {

  /// @class    MethodCloseFiles
  /// @ingroup  Problem
  /// @brief    [\ref Problem] Close any open files

public: // interface

  /// Constructor
  MethodCloseFiles(double seconds_stagger,
		   double seconds_delay,
		   int group_size) throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodCloseFiles);
    
  /// Charm++ PUP::able migration constructor
  MethodCloseFiles (CkMigrateMessage *m)
    : Method (m)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
    TRACEPUP;
    Method::pup(p);
    p | seconds_stagger_;
    p | seconds_delay_;
    p | group_size_;
  };
  
public: // virtual methods

  /// Apply the method 
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "close_files"; }

private: // functions

  void throttle_stagger_();
  void throttle_delay_();
  
private: // attributes

  double seconds_stagger_;
  double seconds_delay_;
  int group_size_;
  
public: // static attributes

  static CmiNodeLock node_lock;

};

#endif /* PROBLEM_METHOD_CLOSE_FILES_HPP */

