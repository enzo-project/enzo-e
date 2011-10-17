// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcessSerial.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Aug 11 12:34:18 PDT 2010
/// @brief    [\ref Parallel] Interface for the GroupProcessSerial class

#ifndef PARALLEL_GROUP_PROCESS_SERIAL_HPP
#define PARALLEL_GROUP_PROCESS_SERIAL_HPP

class Reduce;

class GroupProcessSerial : public GroupProcess {

  /// @class    GroupProcessSerial
  /// @ingroup  Parallel
  /// @brief    [\ref Parallel] Serial implementation of GroupProcess

public: // interface

  /// Initialize an GroupProcessSerial object
  GroupProcessSerial() throw()
  {};

public: // interface (Group)

  /// Synchronize between all compute elements in the Group
  void barrier()  throw()
  {};

  /// Synchronize between two compute elements in the Group
  void sync (int rank, int tag=0) throw()
  {};

  //--------------------------------------------------

  /// Initiate sending an array
  void * send_begin (int    rank_dest, 
		     void * buffer, 
		     int    size, 
		     int    tag=0) throw();

  /// Clean up after sending an array
  void send_end(void * handle) throw()
  {};

  /// Initiate receiving an array
  void * recv_begin
  (int rank_source, void * buffer, int size, int tag=0) throw()
  {
    buffer = buffer_[(long int)tag];
    return (void *) buffer;
  };

  /// Clean up after receiving an array
  void recv_end(void * handle) throw();

  /// Send and receive data between two processes
  void send_recv
  (int rank, void * buffer, int size, int tag=0) throw();

  /// Complete sending or receiving an array
  void wait(void * handle) throw()
  {}

  /// Test completeness of sending or receiving an array
  bool test (void * handle) throw()
  {return true;}

  //--------------------------------------------------

  /// Add an array to a list of arrays to send in bulk
  void bulk_send_add(int rank_dest, void * buffer, int size, int tag=0) throw()
  {};

  /// Initiate a bulk send of multiple arrays
  void * bulk_send() throw()
  {return 0;};

  /// Add an array to a list of arrays to receive in bulk
  void bulk_recv_add(int rank_source, void * buffer, int size, int tag=0) throw()
  {};

  /// Initiate a bulk receive of multiple arrays
  void * bulk_recv() throw()
  {return 0;};

  /// Complete a bulk send or receive of multiple arrays
  void bulk_wait(void * handle) throw()
  {};

  /// Create a Reduce object for this ProcessGroup
  Reduce * create_reduce () throw ();

private: // attributes

  std::map<long int,void *> buffer_;  // Mapping from tag to buffer

};

#endif /* PARALLEL_GROUP_PROCESS_SERIAL_HPP */

