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
    p | offset_;
    p | length_;
  }

  /// Reserve space for a new scalar
  int new_value (std::string name, int n=1)
  {
    int index = name_.size();
    name_.push_back(name);
    int offset = (offset_.size() == 0) ?
      0 : offset_[index-1]+length_[index-1];
    offset_.push_back(offset);
    length_.push_back(n);
    return index;
  }

  /// Return the *last* index of the named scalar
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

  /// Return the offset of items of the given scalar data
  int offset (int index) const { return offset_[index]; }

  /// Return the number of items of the given scalar data (default 1)
  int length (int index) const { return length_[index]; }

  /// Return the number of variables currently stored (including
  /// scalar array elements)
  int size() const {
    int index = name_.size() - 1;
    return (index < 0) ? 0 : (offset_[index] + length_[index]);
  }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Vector of scalar names
  std::vector <std::string> name_;

  /// Vector of offsets of first value
  std::vector<int> offset_;

  /// Vector of scalar data lengths
  std::vector <int> length_;
};

#endif /* DATA_SCALARDESCR_HPP */
