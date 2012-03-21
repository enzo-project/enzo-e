// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputData.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2011-09-19
/// @brief    [\ref Io] Declaration of the OutputData class

#ifndef IO_OUTPUT_DATA_HPP
#define IO_OUTPUT_DATA_HPP

class Factory;
class Hierarchy;
class Patch;
class FieldDescr;

class OutputData : public Output {

  /// @class    OutputData
  /// @ingroup  Io
  /// @brief    [\ref Io] define interface for data I/O

public: // functions

  /// Create an uninitialized OutputData object
  OutputData(const Factory * factory) throw();

  /// Create an uninitialized OutputData object
  ~OutputData() throw()
  {};

public: // virtual functions

  /// Open (or create) a file for IO
  virtual void open () throw();

  /// Close file for IO
  virtual void close () throw();

  /// Finalize output
  virtual void finalize () throw ();

  /// Write hierarchy data to disk
  virtual void write_hierarchy
  ( const Hierarchy * hierarchy,
    const FieldDescr * field_descr ) throw();

  /// Write patch data to disk
  virtual void write_patch
  ( const Patch * patch,
    const FieldDescr * field_descr,
    int ixp0=0, int iyp0=0, int izp0=0) throw();

#ifdef CONFIG_USE_CHARM
  /// Cleanup after writing blocks in a patch
  virtual void end_write_patch () throw();
#endif

  /// Write block data to disk
  virtual void write_block
  ( const Block * block,
    const FieldDescr * field_descr,
    int ixp0=0, int iyp0=0, int izp0=0) throw();

  /// Write local field to disk
  virtual void write_field
  ( const FieldBlock * field_block,
    const FieldDescr * field_descr,
    int field_index) throw();

protected:

};

#endif /* IO_OUTPUT_DATA_HPP */
