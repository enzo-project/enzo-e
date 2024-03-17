// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodCheck.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-02-12
/// @brief    [\ref Enzo] Implementation of EnzoMethodCheck

#ifndef ENZO_ENZO_METHOD_CHECK_HPP
#define ENZO_ENZO_METHOD_CHECK_HPP

class EnzoMethodCheck : public Method {

  /// @class    EnzoMethodCheck
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Write Enzo-E Checkpoint files

public: // interface

  /// Create a new EnzoMethodCheck object
  EnzoMethodCheck
  (int num_files, std::string ordering,
   std::vector<std::string> directory,
   int monitor_iter,
   bool include_ghosts);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodCheck);

  /// Charm++ PUP::able migration constructor
  EnzoMethodCheck (CkMigrateMessage *m)
    : Method (m)
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

public: // virtual methods

  /// Apply the method to advance a block one timestep
  virtual void compute( Block * block) throw();

  virtual std::string name () throw ()
  { return "check"; }

protected: // methods

  DataMsg * create_data_msg_ (Block * block);

protected: // attributes

  /// Number of files to write EnzoBlock data to
  int num_files_;

  /// Ordering method for EnzoBlocks; default "order_morton"
  std::string ordering_;

  /// Disk directory for writing checkpoint files
  std::vector<std::string> directory_;

  /// Whether to include ghosts zones
  bool include_ghosts_;
};

#endif /* ENZO_ENZO_METHOD_CHECK_HPP */
