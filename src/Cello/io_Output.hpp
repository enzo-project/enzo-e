// See LICENSE_CELLO file for license and copyright information

/// @file     io_Output.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2011-03-14
/// @brief    [\ref Io] Declaration of the Output class

#ifndef IO_OUTPUT_HPP
#define IO_OUTPUT_HPP

// #define DEBUG_OUTPUT

class Factory;
class Hierarchy;
class ItIndex;
class ItParticle;
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
  Output (CkMigrateMessage *m)
    : PUP::able(m),
      file_(nullptr),           // Initialization deferred
      schedule_(nullptr),
      sync_write_(1),     // default process-per-stride
      index_(0),
      cycle_(0),
      count_(0),
      time_(0),
      file_name_(""),     // set_filename()
      file_args_(),       // set_filename()
      dir_name_(""),     // set_dirname()
      dir_args_(),
      io_block_(nullptr),
      io_field_data_(nullptr),
      io_particle_data_(nullptr),
      it_field_index_(nullptr),        // set_it_index_field()
      it_particle_index_(nullptr),        // set_it_index_particle()
      stride_write_(1),// default one file per process
      stride_wait_(0) // default no synchronization of writes
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Set file name
  void set_filename (std::string file_name,
		     std::vector<std::string> file_args) throw()
  {
    file_name_ = file_name;
    file_args_ = file_args;
  }

  /// Set dir name
  void set_dir (std::string dir_name,
                std::vector<std::string> dir_args) throw()
  {
    dir_name_ = dir_name;
    dir_args_ = dir_args;
  }

  /// Return the IoBlock object
  IoBlock * io_block () const throw()
  { return io_block_; }

  /// Return the IoFieldData object
  IoFieldData * io_field_data () const throw()
  { return io_field_data_; }

  /// Return the IoParticleData object
  IoParticleData * io_particle_data () const throw()
  { return io_particle_data_; }

  /// Set field iterator
  void set_it_field_index (ItIndex * it_index) throw()
  { it_field_index_ = it_index; }

  /// Set particle iterator
  void set_it_particle_index (ItIndex * it_index) throw()
  { it_particle_index_ = it_index; }
  
  /// Return the File object pointer
  File * file() throw() 
  { return file_; }

  /// Return the Schedule object pointer
  Schedule * schedule() throw() 
  { return schedule_; }

  /// Set schedule
  void set_schedule (Schedule * schedule) throw();

  /// Return whether output is scheduled for this cycle
  bool is_scheduled (int cycle, double time) throw();

  void set_stride_write (int stride) throw () 
  {
    stride_write_ = stride; 
    sync_write_.set_stop(stride_write_);
  }

  int stride_write () const throw () 
  { return stride_write_; }

  int stride_wait () const throw () 
  { return stride_wait_; }

  /// Return whether this process is a writer
  bool is_writer () const throw () 
  { return (CkMyPe() == process_writer()); }

  /// Return the process id of the writer for this process id
  int process_writer() const throw()
  {
    const int ip=CkMyPe();
    return ip - (ip % stride_write_);
  }

  /// Return the updated timestep if time + dt goes past a scheduled output
  double update_timestep (double time, double dt) const throw ();

  /// Write metadata to the file
  void write_meta ( Io * io ) throw ()
  { write_meta_ (meta_type_file, io); }

  /// Write metadata to the current group in the file
  void write_meta_group ( Io * io ) throw ()
  { write_meta_ (meta_type_group, io); }

  /// Accessor function for Charm synchronization of writers
  Sync * sync_write () { return & sync_write_; }

  /// Return the index id in the containing Problem
  int index() const throw() { return index_; }

  /// Advance to next scheduled output
  void next() throw();

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
  {
#ifdef DEBUG_OUTPUT
    CkPrintf ("%d TRACE_OUTPUT Output::write_simulation()\n",CkMyPe());
#endif    
    write_simulation_(simulation);
  }

  /// Write Hierarchy data to disk
  virtual void write_hierarchy ( const Hierarchy * hierarchy ) throw()
  {
#ifdef DEBUG_OUTPUT
    CkPrintf ("%d TRACE_OUTPUT Output::write_hierarchy()\n",CkMyPe());
#endif    
    write_hierarchy_(hierarchy);
  }

  /// Write local block data to disk
  virtual void write_block ( const Block * block ) throw()
  {
#ifdef DEBUG_OUTPUT
    CkPrintf ("%d TRACE_OUTPUT Output::write_block()\n",CkMyPe());
#endif    
    write_block_(block);
  }

  /// Write local field to disk
  virtual void write_field_data ( const FieldData * field_data,
				  int index_field) throw() = 0;

  /// Write local particle data to disk
  virtual void write_particle_data ( const ParticleData * particle_data,
				     int particle_index) throw() = 0;

  /// Prepare local array with data to be sent to remote chare for processing
  virtual void prepare_remote (int * n, char ** buffer) throw()
  {}

  /// Accumulate and write data sent from a remote processes
  virtual void update_remote  ( int n, char * buffer) throw()
  {}

  /// Free local array if allocated; NOP if not
  virtual void cleanup_remote (int * n, char ** buffer) throw()
  {}

protected:

  /// Return the name for the format and given arguments
  std::string expand_name_
  (const std::string * file_name,
   const std::vector<std::string> * file_args) const throw();

  /// Return the path for this file group output.  Creates
  /// the subdirectories if they don't exist
  std::string directory () const
  {
    std::string dir = ".";
    std::string name_dir = expand_name_(&dir_name_,&dir_args_);

    // Create subdirectory if any
    if (name_dir != "") {
      dir = name_dir;
      boost::filesystem::path directory(name_dir);

      if (! boost::filesystem::is_directory(directory)) {

        boost::filesystem::create_directory(directory);

	ASSERT1 ("Output::directory()",
		 "Error creating directory %s",
		 name_dir.c_str(),
                 boost::filesystem::is_directory(directory));
      }
    }

    return dir;
  }

private:

  /// "Loop" over writing the Hierarchy in the Simulation
  void write_simulation_ (const Simulation * simulation ) throw();

  /// Loop over writing Blocks in the Hierarchy
  void write_hierarchy_ ( const Hierarchy * hierarchy ) throw();

  /// Loop over writing Field data in the Block
  void write_block_ ( const Block * block ) throw();

  /// Implementation of write_meta() and write_meta_group()
  void write_meta_ ( meta_type type, Io * io ) throw();

protected: // attributes

  /// File object for output
  File * file_;

  /// Scheduler for this output
  Schedule * schedule_;

  /// Sync for waiting for writers
  Sync sync_write_;

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

  /// Name of the directory, including format arguments 
  std::string dir_name_;

  /// Format strings for dir_name_, if any ("cycle", "time", etc.)
  std::vector<std::string> dir_args_;

  /// I/O Block data accessor
  IoBlock * io_block_;

  /// I/O FieldData data accessor
  IoFieldData * io_field_data_;

  /// I/O ParticleData data accessor
  IoParticleData * io_particle_data_;

  /// Iterator over field indices
  ItIndex * it_field_index_;

  /// Iterator over particle type indices
  ItIndex * it_particle_index_;

  /// Only processes with id's divisible by stride_write_ writes
  /// (1: all processes write; 2: 0,2,4,... write; np: root process writes)
  int stride_write_;
  
  int stride_wait_;

};

#endif /* IO_OUTPUT_HPP */
