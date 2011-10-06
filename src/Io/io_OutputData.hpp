// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputData.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2011-09-19
/// @brief    [\ref Io] Declaration of the OutputData class

#ifndef IO_OUTPUT_DATA_HPP
#define IO_OUTPUT_DATA_HPP

class Hierarchy;
class Patch;
class Schedule;

class ItField;

class OutputData : public Output {

  /// @class    OutputData
  /// @ingroup  Io
  /// @brief    [\ref Io] define interface for data I/O

public: // functions

  /// Create an uninitialized OutputData object
  OutputData(Simulation * simulation) throw();

  /// OutputData destructor
  virtual ~OutputData() throw();

public: // virtual functions

  /// Initialize next output
  virtual void init () throw();

  /// Open (or create) a file for IO
  virtual void open () throw();

  /// Close file for IO
  virtual void close () throw();

  /// Finalize output
  virtual void finalize () throw ();


  /// Write hierarchy data to disk
  virtual void write_hierarchy
  ( const FieldDescr * field_descr,
    Hierarchy * hierarchy) throw();

  /// Write patch data to disk
  virtual void write_patch
  ( const FieldDescr * field_descr,
    Patch * patch,
    int ixp0=0, int iyp0=0, int izp0=0) throw();

  /// Write block data to disk
  virtual void write_block
  ( const FieldDescr * field_descr,
    Block * block,
    int ixp0=0, int iyp0=0, int izp0=0) throw();

  /// Prepare local array with data to be sent to remote chare for processing
  virtual void prepare_remote (int * n, char ** buffer) throw();

  /// Accumulate and write data sent from a remote processes
  virtual void update_remote  ( int n, char * buffer) throw();

  /// Free local array if allocated; NOP if not
  virtual void cleanup_remote (int * n, char ** buffer) throw();


};

#endif /* IO_OUTPUT_DATA_HPP */
