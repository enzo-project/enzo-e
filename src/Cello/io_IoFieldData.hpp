// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoFieldData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoFieldData class

#ifndef IO_IO_FIELD_DATA_HPP
#define IO_IO_FIELD_DATA_HPP

class FieldData;

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

  /// CHARM++ PUP::able declaration
  PUPable_decl(IoFieldData);

  /// CHARM++ migration constructor
  IoFieldData(CkMigrateMessage *m) : Io(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Set FieldIndex
  void set_field_index (int field_index) throw()
  { field_index_ = field_index;};

  /// Set FieldData
  void set_field_data (FieldData * field_data) throw()
  { field_data_ = field_data;};


  /// Return the ith data item associated with the object
  virtual void field_array 
  (void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0,
   int * nx=0,  int * ny=0,  int * nz=0) throw();

  /// Return the ith metadata item associated with the object
  virtual void meta_value 
  (int index, 
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Copy the values to the object
  virtual void save_to (void *); 

protected: // functions

  /// Current FieldData
  FieldData * field_data_;

  /// Index of the field in the FieldData
  int field_index_;


private: // attributes


};

#endif /* IO_IO_FIELD_DATA_HPP */

