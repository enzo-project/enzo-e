// See LICENSE_CELLO file for license and copyright information

/// @file     io_Output.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @todo     Extract Schedule object for scheduling output
/// @brief    [\ref Io] Declaration of the Output class

#ifndef IO_OUTPUT_HPP
#define IO_OUTPUT_HPP

class Hierarchy;
class Patch;
class Schedule;

class ItField;

class Output {

  /// @class    Output
  /// @ingroup  Io
  /// @brief [\ref Io] define interface for various types of IO for
  /// Simulations

public: // functions

  /// Create an uninitialized Output object
  Output(Simulation * simulation ) throw();

  /// Set file name
  void set_filename (std::string filename,
		     std::vector<std::string> fileargs) throw();

  /// Set field iterator
  void set_it_field (ItField * it_field) throw()
  { it_field_ = it_field; }
  

  /// Return the File object pointer
  File * file() throw() 
  { return file_; };

  /// Return the Schedule object pointer
  Schedule * schedule() throw() 
  { return schedule_; };

  /// Return the filename for the file format and given arguments
  std::string expand_file_name () const throw();

  int process_stride () const throw () 
  { return process_stride_; };

  /// Return whether output is scheduled for this cycle
  bool is_scheduled (int cycle, double time);

  /// Return whether this process is a writer
  bool is_writer () const throw () 
  { return (process_ % process_stride_ == 0); };

  double update_timestep (double time, double dt) const;

#ifdef CONFIG_USE_CHARM

  int counter() 
  { 
    if (count_reduce_ >= process_stride_) {count_reduce_ = 0;}
    count_reduce_++; 
    PARALLEL_PRINTF ("Output:counter(%d)\n",count_reduce_);
    return (count_reduce_);
  }

#endif

public: // virtual functions

  /// prepare for accumulating block data
  virtual void init () throw() = 0;

  /// Open (or create) a file for IO
  virtual void open () throw() = 0;

  /// Close file for IO
  virtual void close () throw() = 0;

  /// Write hierarchy data to disk
  virtual void write 
  ( const FieldDescr * field_descr,
    Hierarchy * hierarchy ) throw() = 0;

  /// Write patch data to disk
  virtual void write 
  ( const FieldDescr * field_descr,
    Patch * patch,
    int ix0=0, int iy0=0, int iz0=0) throw() = 0;

  /// Write block data to disk
  virtual void write 
  ( const FieldDescr * field_descr,
    Block * block,
    int ix0=0, int iy0=0, int iz0=0) throw() = 0;

protected: // attributes

  /// Parent simulation (for accessing Hierarchy, etc.)
  Simulation * simulation_;

  /// File object for output
  File * file_;

  /// Scheduler for this output
  Schedule * schedule_;

  /// Only processes with id's divisible by process_stride_ writes
  /// (1: all processes write; 2: 0,2,4,... write; np: root process writes)
  int process_stride_;

  /// ID of this process
  int process_;

#ifdef CONFIG_USE_CHARM

  /// counter for reduction of data from non-writers to writers
  int count_reduce_;

#endif

  /// Simulation cycle for next IO
  int cycle_;

  /// Simulation time for next IO
  double time_;

  /// Name of the file to write, including format arguments
  std::string file_name_;

  /// Format strings for file name, if any ("cycle", "time", etc.)
  std::vector<std::string> file_args_;

  /// Iterator over field id's
  ItField * it_field_;


};

#endif /* IO_OUTPUT_HPP */
