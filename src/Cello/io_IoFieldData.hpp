// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoFieldData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoFieldData class

#ifndef IO_IO_FIELD_DATA_HPP
#define IO_IO_FIELD_DATA_HPP

class FieldData;
class FieldDescr;

class IoFieldData : public Io {

  /// @class    IoFieldData
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between FieldData and Output classes

public: // interface

  /// Constructor
  IoFieldData() throw();

  /// Destructor
  virtual ~IoFieldData() throw()
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Set FieldIndex
  void set_field_index (int field_index) throw()
  { field_index_ = field_index;};

  /// Set FieldData
  void set_field_data (FieldData * field_data) throw()
  { field_data_ = field_data;};

  /// Set FieldDescr
  void set_field_descr (FieldDescr * field_descr) throw()
  { field_descr_ = field_descr; };


#include "_io_Io_common.hpp"

  
protected: // functions

  /// FieldDescr for the FieldData
  FieldDescr * field_descr_;

  /// Current FieldData
  FieldData * field_data_;

  /// Index of the field in the FieldData
  int field_index_;


private: // attributes


};

#endif /* IO_IO_FIELD_DATA_HPP */

