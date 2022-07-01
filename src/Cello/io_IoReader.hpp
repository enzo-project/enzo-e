// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoReader.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-02-01
/// @brief    [\ref Io] Declaration of the IoReader class

#ifndef IO_IO_READER_HPP
#define IO_IO_READER_HPP

class IoReader : public CBase_IoReader {

  /// @class    IoReader
  /// @ingroup  Io
  /// @brief    [\ref Io] 

public: // interface

  /// Constructor
  IoReader() throw()
  : CBase_IoReader()
  { }

  /// CHARM++ migration constructor
  IoReader (CkMigrateMessage *m) :
    CBase_IoReader(m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  {
    TRACEPUP;
    CBase_IoReader::pup(p); 
  }

private: // functions

private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* IO_IO_READER_HPP */

