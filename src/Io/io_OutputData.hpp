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

#ifdef CONFIG_USE_CHARM

  /// Prepare for accumulating block data
  virtual void init (const Hierarchy * hierarchy, int cycle, double time) throw();

  /// Accumulate block-local data
  virtual void block (const Block * block) throw();

#endif

  /// Write hierarchy data to disk
  virtual void write 
  ( const FieldDescr * field_descr,
    Hierarchy * hierarchy, 
    int cycle, double time,
    bool root_call=true) throw();

  /// Write patch data to disk
  virtual void write 
  ( const FieldDescr * field_descr,
    Patch * patch,
    int cycle, double time, 
    bool root_call=true, int ix0=0, int iy0=0, int iz0=0) throw();

  /// Write block data to disk
  virtual void write 
  ( const FieldDescr * field_descr,
    Block * block,
    int cycle, double time, 
    bool root_call=true, int ix0=0, int iy0=0, int iz0=0) throw();

};

#endif /* IO_OUTPUT_DATA_HPP */
