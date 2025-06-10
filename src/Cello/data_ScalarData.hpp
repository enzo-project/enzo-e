// See LICENSE_CELLO file for license and copyright information

/// @file     data_ScalarData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2018-07-16
/// @brief    [\ref Data] Declaration of the ScalarData class

#ifndef DATA_SCALARDATA_HPP
#define DATA_SCALARDATA_HPP

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

  /// Return the specified scalar value
  T * value (const ScalarDescr * scalar_descr, int index)
  { if (index >= int(value_.size())) allocate(scalar_descr);
    return (index >= 0) ? & value_[scalar_descr->offset(index)] : NULL;
  }

  const T * value (const ScalarDescr * scalar_descr, int index) const 
  { return (index >= 0) ? & value_.at(scalar_descr->offset(index)) : NULL; }

  /// Allocate scalars
  void allocate(const ScalarDescr * scalar_descr)
  { value_.resize(scalar_descr->size()); }

  size_t size() const
  { return value_.size(); }
  
  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const
  {
    int size = 0;
    SIZE_VECTOR_TYPE(size,T,value_);
    return size;
  };

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const
  {
    char * pc = buffer;
    SAVE_VECTOR_TYPE(pc,T,value_);
    return pc;
  }

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer)
  {
    char * pc = buffer;
    LOAD_VECTOR_TYPE(pc,T,value_);
    return pc;
  }

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Vector of values
  std::vector <T> value_;
};

#endif /* DATA_SCALARDATA_HPP */

