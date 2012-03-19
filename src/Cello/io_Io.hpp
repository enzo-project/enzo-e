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
  Io(size_t meta_count, size_t data_count) throw();

  /// Destructor
  virtual ~Io () throw ()
  {}


#include "_io_Io_common.hpp"


  /// Return number of metadata items associated with the associated class
  size_t meta_count() const throw()
  { return meta_count_; }

  /// Return number of data items associated with the associated class
  size_t data_count() const throw()
  { return data_count_; }

  
protected: // attributes

  /// Number of metadata items
  size_t meta_count_;

  /// Name of the metadata items
  std::vector <std::string> meta_name_;

  /// Number of data items
  size_t data_count_;

  /// Name of the data items
  std::vector <std::string> data_name_;


};

#endif /* IO_IO_HPP */

