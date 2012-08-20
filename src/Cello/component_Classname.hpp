// See LICENSE_CELLO file for license and copyright information

/// @file     component_Classname.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Component] Declaration of the Classname class
///

#ifndef COMPONENT_CLASSNAME_HPP
#define COMPONENT_CLASSNAME_HPP

class Classname {

  /// @class    Classname
  /// @ingroup  Component
  /// @brief    [\ref Component] 

public: // interface

  /// Constructor
  Classname() throw();

  /// Destructor
  ~Classname() throw();

  /// Copy constructor
  Classname(const Classname & classname) throw();

  /// Assignment operator
  Classname & operator= (const Classname & classname) throw();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change
  }
#endif
  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* COMPONENT_CLASSNAME_HPP */

