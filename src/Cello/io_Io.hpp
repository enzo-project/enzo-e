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
  Io(size_t data_count = 0) throw();

  /// Destructor
  virtual ~Io () throw ()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    p | meta_name_;
    p | data_count_;

  }

#include "_io_Io_common.hpp"

  /// Return number of metadata items associated with the associated class
  size_t meta_count() const throw()
  { return meta_name_.size(); }

  /// Return number of data items associated with the associated class
  size_t data_count() const throw()
  { return data_count_; }

  
protected: // attributes

  /// Name of the metadata items
  std::vector <std::string> meta_name_;

  /// Name of the data items
  size_t data_count_;


};

#endif /* IO_IO_HPP */

