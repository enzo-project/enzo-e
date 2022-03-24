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
  : CBase_IoEnzoReader(),
    name_dir_(),
    name_file_(),
    stream_block_list_(),
    file_(nullptr)
  { }

  IoEnzoReader(std::string name_dir,
               std::string name_file) throw();

  /// CHARM++ migration constructor
  IoEnzoReader(CkMigrateMessage *m) : CBase_IoEnzoReader(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; CBase_IoEnzoReader::pup(p); }

  
protected: // functions

  FileHdf5 * file_open_(std::string name_dir, std::string name_file);
  std::ifstream open_block_list_(std::string name_dir, std::string name_file);
  void file_read_hierarchy_();
  void file_read_block_(EnzoMsgCheck * msg_check, std::string file_name);
  void read_meta_ ( FileHdf5 * file, Io * io, std::string type_meta );
  bool read_block_list_(std::string & block_name);
  void close_block_list_();
  void read_dataset_
  (File * file, char * buffer, int type_data,
   int nx, int ny, int nz,
   int m4[4]);

  template <class T>
  void copy_buffer_to_particle_attribute_
  (T * buffer, Particle particle, int it, int ia, int np);


private: // attributes

  // NOTE: change pup() function whenever attributes change
  std::string name_dir_;
  std::string name_file_;

  std::ifstream stream_block_list_;

  FileHdf5 * file_;
};

#endif /* ENZO_IO_ENZO_READER_HPP */

