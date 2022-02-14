// See LICENSE_CELLO file for license and copyright information

/// @file     io_Io.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Abstract base class for providing access to class data for Output classes

#ifndef IO_IO_HPP
#define IO_IO_HPP

class Io : public PUP::able {

  /// @class    Io
  /// @ingroup  Io
  /// @brief    [\ref Io] Declaration of the Io base class

public: // interface

  /// Constructor
  Io() throw();

  /// Destructor
  virtual ~Io () throw ()
  {}

  /// CHARM++ function for determining the specific class in the class hierarchy
  PUPable_decl(Io);

  Io (CkMigrateMessage *m) : PUP::able(m) 
  { TRACE("Io::Io(CkMigrateMessage*)"); }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    PUP::able::pup(p);
    TRACEPUP;

    p | meta_name_;

  }

  /// PACKING / UNPACKING
  
  /// Return the number of bytes required to serialize the data object
  virtual int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  virtual char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  virtual char * load_data (char * buffer);

#include "_io_Io_common.hpp"

  /// Return number of metadata items associated with the associated class
  size_t meta_count() const throw()
  { return meta_name_.size(); }

protected: // attributes

  /// Name of the metadata items
  std::vector <std::string> meta_name_;

};

#endif /* IO_IO_HPP */

