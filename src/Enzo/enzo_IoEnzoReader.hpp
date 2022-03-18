// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_IoEnzoReader.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-02-01
/// @brief    [\ref Io] Declaration of the IoEnzoReader class

#ifndef ENZO_IO_ENZO_READER_HPP
#define ENZO_IO_ENZO_READER_HPP

class IoEnzoReader : public CBase_IoEnzoReader {

  /// @class    IoEnzoReader
  /// @ingroup  Io
  /// @brief    [\ref Io] 

public: // interface

  /// Constructors
  IoEnzoReader() throw()
    : CBase_IoEnzoReader()
  { CkPrintf ("%d IoEnzoReader()\n",CkMyPe()); }

  IoEnzoReader(std::string file_name) throw();

  /// CHARM++ migration constructor
  IoEnzoReader(CkMigrateMessage *m) : CBase_IoEnzoReader(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; CBase_IoEnzoReader::pup(p); }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change
  std::string file_name_;

};

#endif /* ENZO_IO_ENZO_READER_HPP */

