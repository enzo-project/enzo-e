// See LICENSE_CELLO file for license and copyright information

/// @file     control_Control.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-11-03 23:12:52
/// @brief    [\ref Control] Declaration of the Control class
///

#ifndef CONTROL_CONTROL_HPP
#define CONTROL_CONTROL_HPP

// class Control : public PUP::able {

//   /// @class    Control
//   /// @ingroup  Control
//   /// @brief    [\ref Control] 

// public: // interface

//   /// empty constructor for charm++ pup()
//   Control() throw() {};

//   /// Constructor
//   Control(std::string curr_name,
// 	  std::string next_name,
// 	  int phase, std::string sync) throw() 
//     : phase_(phase),
//       sync_(sync),
//       curr_name_(curr_name),
//       next_name_(next_name)
//   {}

//   /// Copy constructor
//   Control(const Control & Control) throw()
//   {}

//   /// CHARM++ PUP::able declaration
//   PUPable_decl(Control);

//   /// CHARM++ migration constructor for PUP::able

//   Control (CkMigrateMessage *m) : PUP::able(m) {}

//   /// CHARM++ Pack / Unpack function
//   void pup (PUP::er &p)
//   { 
//     // NOTE: change this function whenever attributes change
//     PUP::able::pup(p);
//     TRACEPUP;

//     p | phase_;
//     p | sync_;
//     p | curr_name_;
//     p | next_name_;
//     //    p | next_func_;
//     // WARNING("Control::pup()",
//     // 	    "skipping next_func_");
//   }

//   /// Return the phase of the calculation associated with this Control object
//   int phase() const
//   { return phase_; }

//   /// Return the type of synchronization associated with this Control object
//   std::string sync() const
//   { return sync_; }

//   /// Return the name of the function that created this Control object
//   std::string curr_name() const
//   { return curr_name_; }

//   /// Return the name of the function that this Control object will call next
//   std::string next_name() const
//   { return next_name_; }

//   void print () const
//   {
//     printf ("%s:%d DEBUG  %s calling %s phase %d using %s synchronization\n",
// 	    __FILE__,__LINE__,
// 	    curr_name_.c_str(),
// 	    next_name_.c_str(),
// 	    phase_,
// 	    sync_.c_str());
//   }
  
// private: // attributes

//   int phase_;
//   std::string sync_;
//   std::string curr_name_;
//   std::string next_name_;
//   //  void (Block::*next_func_)();

// };

#endif /* CONTROL_CONTROL_HPP */

