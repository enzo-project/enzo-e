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
  IoEnzoReader();
  
  /// CHARM++ migration constructor
  IoEnzoReader(CkMigrateMessage *m) : CBase_IoEnzoReader(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP;
    CBase_IoEnzoReader::pup(p);
    p | name_dir_;
    p | name_file_;
    p | max_level_;
    p | level_;
    //    p | stream_block_list_;
    //    p | file_;
    p | sync_blocks_;
    //    p | io_msg_check_;
  }

  /// Send data to existing root blocks
  void p_init_root
  (std::string name_dir, std::string name_file, int max_level);

  /// Create blocks in the given level
  void p_create_level(int level);

  /// Receive acknowledgement that a block was created

  void p_block_created();

  /// Initialize existing blocks in the given level
  void p_init_level(int level);

  /// Received acknowledgement that the block is done
  void p_block_ready();

protected: // functions

  /// update synchronization given that the given block is done
  void block_ready_();
  void block_created_();

  void file_open_block_list_(std::string name_dir, std::string name_file);
  void file_read_block_(EnzoMsgCheck * msg_check, std::string file_name);
  void file_read_block_fields_(DataMsg * data_msg, int nx, int ny, int nz);
  void file_read_block_particles_(DataMsg * data_msg);
  bool read_block_list_(std::string & block_name, int & level);
  void file_close_block_list_();

  std::ifstream stream_open_blocks_(std::string name_dir, std::string name_file);
  void file_read_hierarchy_();
  void read_meta_ ( FileHdf5 * file, Io * io, std::string type_meta );
  void file_read_dataset_
  (char * buffer, int type_data,
   int nx, int ny, int nz,
   int m4[4]);

  template <class T>
  void copy_buffer_to_particle_attribute_
  (T * buffer, Particle particle, int it, int ia, int np);

private: // attributes

  // NOTE: change pup() function whenever attributes change
  
  std::string name_dir_;
  std::string name_file_;

  int max_level_;
  std::ifstream stream_block_list_;

  FileHdf5 * file_;

  Sync sync_blocks_;

  /// List of blocks in the file by level (negative blocks included in
  /// level 0
  std::vector< std::vector<EnzoMsgCheck *> > io_msg_check_;

  /// Current level
  int level_;
  /// List of block names in the file
  std::vector<std::string> block_name_list_;

  /// List of levels for the in block_name_list_
  std::vector<int>         block_level_list_;

  /// Count of blocks in each level_
  std::vector<int> blocks_in_level_;
};

#endif /* ENZO_IO_ENZO_READER_HPP */

