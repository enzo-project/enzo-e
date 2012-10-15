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
class Config;

class OutputData : public Output {

  /// @class    OutputData
  /// @ingroup  Io
  /// @brief    [\ref Io] define interface for data I/O

public: // functions

  /// Empty constructor for Charm++ pup()
  OutputData() throw() {}

  /// Create an uninitialized OutputData object
  OutputData(int index,
	     const Factory * factory,
	     Config * config) throw();

  /// Close the file if it is open
  virtual ~OutputData() throw();

#ifdef CONFIG_USE_CHARM

  /// Charm++ PUP::able declarations
  PUPable_decl(OutputData);

  /// Charm++ PUP::able migration constructor
  OutputData (CkMigrateMessage *m) : Output (m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#endif

public: // virtual functions

  /// Open (or create) a file for IO
  virtual void open () throw();

  /// Close file for IO
  virtual void close () throw();

  /// Whether the file is open or not
  virtual bool is_open () throw();

  /// Finalize output
  virtual void finalize () throw ();

  /// Write hierarchy data to disk
  virtual void write
  ( const Hierarchy * hierarchy,
    const FieldDescr * field_descr ) throw();

  /// Write patch data to disk
  virtual void write
  ( const Patch * patch,
    const FieldDescr * field_descr,
    int ixp0=0, int iyp0=0, int izp0=0) throw();

  /// Write block data to disk
  virtual void write
  ( const Block * block,
    const FieldDescr * field_descr,
    int ixp0=0, int iyp0=0, int izp0=0) throw();

  /// Write local field to disk
  virtual void write
  ( const FieldBlock * field_block,
    const FieldDescr * field_descr,
    int field_index) throw();

protected:

};

#endif /* IO_OUTPUT_DATA_HPP */
