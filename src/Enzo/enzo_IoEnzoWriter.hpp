// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_IoEnzoWriter.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-02-11
/// @brief    [\ref Io] Declaration of the IoEnzoWriter class

#ifndef ENZO_IO_ENZO_WRITER_HPP
#define ENZO_IO_ENZO_WRITER_HPP

class IoEnzoWriter : public CBase_IoWriter {

  /// @class    IoEnzoWriter
  /// @ingroup  Io
  /// @brief    [\ref Io]

public: // interface

  /// Constructor
  IoEnzoWriter() throw();

  /// CHARM++ migration constructor
  IoEnzoWriter(CkMigrateMessage *m) : CBase_IoWriter(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP; CBase_IoWriter::pup(p); }

public: // entry methods

  void p_write()
  { CkPrintf ("DEBUG_IO IoEnzoWriter::p_write()\n"); }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* ENZO_IO_ENZO_WRITER_HPP */

