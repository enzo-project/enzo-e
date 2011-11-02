// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoFieldBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoFieldBlock class

#ifndef IO_IO_FIELD_BLOCK_HPP
#define IO_IO_FIELD_BLOCK_HPP

class FieldBlock;

class IoFieldBlock : public Io {

  /// @class    IoFieldBlock
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between FieldBlock and Output classes

public: // interface

  /// Constructor
  IoFieldBlock() throw();

  /// Destructor
  virtual ~IoFieldBlock() throw()
  {}

  /// Set FieldIndex
  void set_field_index (int field_index) throw()
  { field_index_ = field_index;};

  /// Set FieldBlock
  void set_field_block (const FieldBlock * field_block) throw()
  { field_block_ = field_block;};

  /// Set FieldDescr
  void set_field_descr (const FieldDescr * field_descr) throw()
  { field_descr_ = field_descr; };

  /// Return the ith metadata item associated with the FieldBlock object
  void meta_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  /// Return the ith data item associated with the FieldBlock object
  void data_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  
protected: // functions

  /// FieldDescr for the FieldBlock
  const FieldDescr * field_descr_;

  /// Current FieldBlock
  const FieldBlock * field_block_;

  /// Index of the field in the FieldBlock
  int field_index_;


private: // attributes


};

#endif /* IO_IO_FIELD_BLOCK_HPP */

