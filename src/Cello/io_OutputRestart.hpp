// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputRestart.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2011-09-26
/// @brief    [\ref Io] Declaration for the OutputRestart class

#ifndef IO_OUTPUT_RESTART_HPP
#define IO_OUTPUT_RESTART_HPP

class ItField;

class OutputRestart : public Output {

  /// @class    OutputRestart
  /// @ingroup  Io
  /// @brief [\ref Io] class for reading and writing Restart files

public: // functions

  /// Create an uninitialized OutputRestart object
  OutputRestart(const Factory * factory) throw();

  /// OutputRestart destructor
  virtual ~OutputRestart() throw();

public: // virtual functions

  /// Prepare for accumulating block data
  virtual void init () throw();

  /// Open (or create) a file for IO
  virtual void open () throw();

  /// Close file for IO
  virtual void close () throw();

  /// Finalize output
  virtual void finalize () throw ();

  /// Write an entire simulation to disk
  virtual void write_simulation
  ( 
   Factory    * factory,
   FieldDescr * field_descr,
   Hierarchy  * hierarchy,
   Simulation * simulation
    ) throw();

  /// Write hierarchy-related field data
  virtual void write_hierarchy
  ( const FieldDescr * field_descr,
    Hierarchy * hierarchy) throw();

  /// Write patch-related field data; may be called by write_hierarchy
  virtual void write_patch
  ( const FieldDescr * field_descr,
    Patch * patch,
    int ixp0=0, int iyp0=0, int izp0=0) throw();

  /// Write block-related field data; may be called by write_patch
  virtual void write_block
  ( const FieldDescr * field_descr,
    Block * block,
    int ixp0=0, int iyp0=0, int izp0=0) throw();

  /// Write local field to disk
  virtual void write_field
  ( const FieldDescr * field_descr,
    FieldBlock * field_block, int field_index) throw();

  /// Prepare local array with data to be sent to remote chare for processing
  virtual void prepare_remote (int * n, char ** buffer) throw();

  /// Accumulate and write data sent from a remote processes
  virtual void update_remote  ( int n, char * buffer) throw();

  /// Free local array if allocated; NOP if not
  virtual void cleanup_remote (int * n, char ** buffer) throw();

private: // functions

private: // attributes

};

#endif /* IO_OUTPUT_RESTART_HPP */
