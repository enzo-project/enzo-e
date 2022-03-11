// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_IoEnzoWriter.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-02-11
/// @brief    [\ref Io] Declaration of the IoEnzoWriter class

#ifndef ENZO_IO_ENZO_WRITER_HPP
#define ENZO_IO_ENZO_WRITER_HPP
class IoEnzoWriter : public CBase_IoEnzoWriter {

  /// @class    IoEnzoWriter
  /// @ingroup  Io
  /// @brief    [\ref Io]

public: // interface

  /// Defaut Constructor
  IoEnzoWriter() throw()
  : CBase_IoEnzoWriter(),
    num_files_(0),
    ordering_("")
  {  }

  /// Constructor
  IoEnzoWriter(int num_files, std::string ordering) throw();

  /// CHARM++ migration constructor
  IoEnzoWriter(CkMigrateMessage *m) : CBase_IoEnzoWriter(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP; CBase_IoEnzoWriter::pup(p);

    p | num_files_;
    p | ordering_;
  }

public: // entry methods

  void p_write(int index_file, std::string name_this, std::string name_next,
               Index index_this, Index index_next, int index,
               bool is_last);

  // void r_created(CkReductionMsg *msg);

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change
  int num_files_;
  std::string ordering_;

};

#endif /* ENZO_IO_ENZO_WRITER_HPP */

