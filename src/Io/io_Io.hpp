// See LICENSE_CELLO file for license and copyright information

/// @file     io_Io.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Abstract base class for providing access to class data for Output classes

#ifndef IO_IO_HPP
#define IO_IO_HPP

class Io {

  /// @class    Io
  /// @ingroup  Io
  /// @brief    [\ref Io] Declaration of the Io base class

public: // interface

  /// Constructor
  Io(int meta_count, int data_count) throw();

  /// Destructor
  virtual ~Io () throw ()
  {}

  /// Return the ith metadata item associated with the associated object

  virtual void meta_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw();

  /// Return the ith data item associated with the associated object
  virtual void data_value 
  (int index, 
   void ** buffer, std::string * name, enum scalar_type * type,
   int * n0=0, int * n1=0, int * n2=0, int * n3=0, int * n4=0) throw() = 0;

  /// Return number of metadata items associated with the associated class
  int meta_count() const throw()
  { return meta_count_; }

  /// Return number of data items associated with the associated class
  int data_count() const throw()
  { return data_count_; }

  
protected: // attributes

  /// Number of metadata items
  int meta_count_;

  /// Name of the metadata items
  std::vector <std::string> meta_name_;

  /// Number of data items
  int data_count_;

  /// Name of the data items
  std::vector <std::string> data_name_;


};

#endif /* IO_IO_HPP */

