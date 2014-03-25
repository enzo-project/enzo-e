// See LICENSE_CELLO file for license and copyright information

/// @file     io_Output.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2011-03-14
/// @brief    [\ref Io] Declaration of the Output class

#ifndef IO_OUTPUT_HPP
#define IO_OUTPUT_HPP

class Factory;
class FieldDescr;
class Hierarchy;
class ItField;
class Schedule;
class Simulation;

class Output : public PUP::able 
{

  /// @class    Output
  /// @ingroup  Io
  /// @brief [\ref Io] define interface for various types of IO for
  /// Simulations

public: // functions

  /// Empty constructor for Charm++ pup()
  Output() throw() { }

  /// Create an uninitialized Output object
  Output(int index, const Factory * factory) throw();

  /// Delete an Output object
  virtual ~Output() throw();

  /// Charm++ PUP::able declarations
  PUPable_abstract(Output);

  /// Charm++ PUP::able migration constructor
  Output (CkMigrateMessage *m) : PUP::able(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Set file name
  void set_filename (std::string filename,
		     std::vector<std::string> fileargs) throw()
  { file_name_ = filename;  file_args_ = fileargs;  }

  /// Set field iterator
  void set_it_field (ItField * it_field) throw()
  { it_field_ = it_field; }
  
  /// Return the IoBlock object
  IoBlock * io_block () const throw() { return io_block_; }

  /// Return the IoFieldBlock object
  IoFieldBlock * io_field_block () const throw() { return io_field_block_; }

  /// Return the File object pointer
  File * file() throw() 
  { return file_; };

  /// Return the Schedule object pointer
  Schedule * schedule() throw() 
  { return schedule_; };

  /// Set schedule
  void set_schedule (Schedule * schedule) throw();

  int process_stride () const throw () 
  { return process_stride_; };

  void set_process_stride (int stride) throw () 
  {
    process_stride_ = stride; 
    sync_.set_stop(process_stride_);
  };

  /// Return whether output is scheduled for this cycle
  bool is_scheduled (int cycle, double time);

  /// Return whether this process is a writer
  bool is_writer () const throw () 
  { return (process_ == process_writer()); };

  /// Return the process id of the writer for this process id
  int process_writer() const throw()
  {
    return process_ - (process_ % process_stride_);
  }

  /// Return the updated timestep if time + dt goes past a scheduled output
  double update_timestep (double time, double dt) const throw ();

  /// Write metadata to the file
  void write_meta ( Io * io ) throw ()
  { write_meta_ (meta_type_file, io); }

  /// Write metadata to the current group in the file
  void write_meta_group ( Io * io ) throw ()
  { write_meta_ (meta_type_group, io); }

  /// Accessor function for the CHARM Sync class
  Sync * sync() { return & sync_; };

  /// Return the index id in the containing Problem
  int index() const throw() { return index_; }


public: // virtual functions

  /// Initialize next output
  virtual void init () throw()
  {} ;

  /// Open (or create) a file for IO
  virtual void open () throw() = 0;

  /// Close file for IO
  virtual void close () throw() = 0;

  /// Finalize output
  virtual void finalize () throw ()
  { count_ ++; }

  /// Write Simulation data to disk
  virtual void write_simulation ( const Simulation * simulation ) throw()
  { write_simulation_(simulation); }

  /// Write Hierarchy data to disk
  virtual void write_hierarchy
  ( const Hierarchy * hierarchy, 
    const FieldDescr * field_descr  ) throw()
  { write_hierarchy_(hierarchy,field_descr); }

  /// Write local block data to disk
  virtual void write_block
  ( const CommBlock * block, 
    const FieldDescr * field_descr) throw()
  { write_block_(block,field_descr); }

  /// Write local field to disk
  virtual void write_field_block
  ( const FieldBlock * field_block, 
    const FieldDescr * field_descr,
    int field_index) throw() = 0;

  /// Prepare local array with data to be sent to remote chare for processing
  virtual void prepare_remote (int * n, char ** buffer) throw()
  {};

  /// Accumulate and write data sent from a remote processes
  virtual void update_remote  ( int n, char * buffer) throw()
  {};

  /// Free local array if allocated; NOP if not
  virtual void cleanup_remote (int * n, char ** buffer) throw()
  {};

protected:

  /// "Loop" over writing the Hierarchy in the Simulation
  void write_simulation_ (const Simulation * simulation ) throw();

  /// Loop over writing Blocks in the Hierarchy
  void write_hierarchy_
  ( const Hierarchy * hierarchy, 
    const FieldDescr * field_descr  ) throw();

  /// Loop over writing Field data in the CommBlock
  void write_block_
  ( const CommBlock * block, 
    const FieldDescr * field_descr) throw();

  /// Return the filename for the file format and given arguments
  std::string expand_file_name_
  (const std::string * file_name,
   const std::vector<std::string> * file_args) const throw();

private:

  /// Implementation of write_meta() and write_meta_group()
  void write_meta_ ( meta_type type, Io * io ) throw();


protected: // attributes

  /// File object for output
  File * file_;

  /// Scheduler for this output
  Schedule * schedule_;

  /// ID of this process
  int process_;

  /// Sync for ending output
  Sync sync_;

  /// Index of this Output object in Simulation
  size_t index_;

  /// Simulation cycle for next IO
  int cycle_;

  /// Count of number of times this Output object performed output
  int count_;

  /// Simulation time for next IO
  double time_;

  /// Name of the file to write, including format arguments
  std::string file_name_;

  /// Format strings for file name, if any ("cycle", "time", etc.)
  std::vector<std::string> file_args_;

  /// Iterator over field id's
  ItField * it_field_;

  /// I/O CommBlock data accessor
  IoBlock * io_block_;

  /// I/O FieldBlock data accessor
  IoFieldBlock * io_field_block_;

  /// Only processes with id's divisible by process_stride_ writes
  /// (1: all processes write; 2: 0,2,4,... write; np: root process writes)
  int process_stride_;


};

#endif /* IO_OUTPUT_HPP */
