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
    ordering_(""),
    stream_block_list_(),
    file_(nullptr),
    monitor_iter_(0)
  {  }

  /// Constructor
  IoEnzoWriter(int num_files,
               std::string ordering,
               int monitor_iter) throw();

  /// CHARM++ migration constructor
  IoEnzoWriter(CkMigrateMessage *m) : CBase_IoEnzoWriter(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP; CBase_IoEnzoWriter::pup(p);

    p | num_files_;
    p | ordering_;
    p | monitor_iter_;
  }

public: // entry methods

  void p_write(EnzoMsgCheck *);

  // void r_created(CkReductionMsg *msg);

protected: // functions

  FileHdf5 * file_open_(std::string name_dir, std::string name_file);
  std::ofstream create_block_list_(std::string name_dir, std::string name_file);
  void file_write_hierarchy_();
  void file_write_block_(EnzoMsgCheck * msg_check);
  void write_meta_ ( FileHdf5 * file, Io * io, std::string type_meta );

  void write_block_list_(std::string block_name, int level);
  void close_block_list_();

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Number of files in which to store Block data
  int num_files_;

  /// Block-ordering used for mapping blocks to files
  std::string ordering_;

  std::ofstream stream_block_list_;

  FileHdf5 * file_;

  /// How often to output write status wrt block indices in first
  /// file; 0 for no output
  int monitor_iter_;
};

#endif /* ENZO_IO_ENZO_WRITER_HPP */

