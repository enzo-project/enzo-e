// See LICENSE_CELLO file for license and copyright information

/// @file     data_ScalarData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2018-07-16
/// @brief    [\ref Data] Declaration of the ScalarData class

#ifndef DATA_SCALARDATA_HPP
#define DATA_SCALARDATA_HPP

#ifdef NEW_SYNC
template <class T>
class ScalarData {

  /// @class    ScalarData
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    p | value_;
  }

  /// Allocate scalars
  void allocate(ScalarDescr * scalar_descr)
  { value_.resize(scalar_descr->size()); }
  
  /// Return the specified scalar value
  T & value (ScalarDescr * scalar_descr, int index)
  { if (index >= value_.size()) allocate(scalar_descr);
    return value_[index];
  }

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Vector of values
  std::vector <T> value_;
};
#endif /* NEW_SYNC */

#endif /* DATA_SCALARDATA_HPP */

