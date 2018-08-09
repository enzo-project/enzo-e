// See LICENSE_CELLO file for license and copyright information

/// @file     data_ScalarDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2018-07-16
/// @brief    [\ref Data] Declaration of the ScalarDescr class

#ifndef DATA_SCALARDESCR_HPP
#define DATA_SCALARDESCR_HPP

class ScalarDescr {

  /// @class    ScalarDescr
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p | name_;
  }

  /// Reserve space for a new scalar
  int new_value (std::string name)
  {
    int index = name_.size();
    name_.push_back(name);
    return index;
  }

  /// Return the *first* index of the named scalar
  int index (std::string name) const
  {
    int i=name_.size();
    while (--i>=0) {
      if (name == name_[i]) break;
    }
    return i;
  }

  /// Return the name of the given scalar
  std::string name (int index) const { return name_[index]; }

  /// Return the number of values currently stored
  int size() const { return name_.size(); }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Vector of scalar names
  std::vector <std::string> name_;

};

#endif /* DATA_SCALARDESCR_HPP */
