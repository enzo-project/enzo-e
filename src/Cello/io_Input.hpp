// See LICENSE_CELLO file for license and copyright information

/// @file     io_Input.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2012-03-23
/// @brief    [\ref Io] Declaration of the Input class

#ifndef IO_INPUT_HPP
#define IO_INPUT_HPP

class Factory;
class FieldDescr;
class Hierarchy;
class ItField;
class Patch;
class Simulation;

class Input {

  /// @class    Input
  /// @ingroup  Io
  /// @brief [\ref Io] define interface for various types of IO for
  /// Simulations

public: // functions

  /// Create an uninitialized Input object
  Input(const Factory * factory) throw();

  /// Delete an Input object
  virtual ~Input() throw();

  /// Set file name
  void set_filename (std::string filename,
		     std::vector<std::string> fileargs) throw();

  /// Set field iterator
  void set_it_field (ItField * it_field) throw()
  { it_field_ = it_field; }
  
  /// Return the IoBlock object
  IoBlock * io_block () const throw() { return io_block_; }

  /// Return the IoFieldBlock object
  IoFieldBlock * io_field_block () const throw() { return io_field_block_; }

  /// Return the File object pointer
  File * file() throw() { return file_; };

  int process_stride () const throw () 
  { return process_stride_; };

  void set_process_stride (int stride) throw () 
  {
    process_stride_ = stride; 
#ifdef CONFIG_USE_CHARM
    counter_.set_value(process_stride_);
#endif
  };

  /// Used to set cycle and time; always returns true
  bool is_scheduled (int cycle, double time);

  /// Return whether this process is a reader
  bool is_reader () const throw () 
  { return (process_ == process_reader()); };

  /// Return the process id of the reader for this process id
  int process_reader() const throw()
  {
    return process_ - (process_ % process_stride_);
  }

#ifdef CONFIG_USE_CHARM

  /// Accessor function for the CHARM Counter class
  Counter * counter() { return & counter_; };

  /// Set the index of this input in its simulation
  void set_index_charm(int index_charm) { index_charm_ = index_charm; }

#endif

public: // virtual functions

  /// Initialize next input
  virtual void init () throw()
  {} ;

  /// Open (or create) a file for IO
  virtual void open () throw() = 0;

  /// Whether the file is open or not
  virtual bool is_open () throw() = 0;

  /// Close file for IO
  virtual void close () throw() = 0;

  /// Finalize input
  virtual void finalize () throw ()
  { count_input_ ++; }

  /// Read metadata to the file
  void read_meta ( Io * io ) throw ()
  { read_meta_ (meta_type_file, io); }

  /// Read metadata to the current group in the file
  void read_meta_group ( Io * io ) throw ()
  { read_meta_ (meta_type_group, io); }

public:
  /// Read an entire simulation to disk
  virtual void read_simulation ( Simulation * simulation) throw();

  /// Read local hierarchy data to disk
  virtual void read_hierarchy
  ( Hierarchy * hierarchy, 
    const FieldDescr * field_descr  ) throw();

  /// Read local patch data to disk
  virtual Patch * read_patch
  ( Patch * patch,
    const FieldDescr * field_descr,  
    int ixp0=0, int iyp0=0, int izp0=0) throw();

#ifdef CONFIG_USE_CHARM
  /// Cleanup after writing blocks in a patch
  virtual void end_read_patch () throw()
  { }
#endif

  /// Read local block data to disk
  virtual Block * read_block
  ( Block * block, 
    const FieldDescr * field_descr,  
    int ixp0=0, int iyp0=0, int izp0=0) throw();

  /// Read local field to disk
  virtual void read_field
  ( FieldBlock * field_block, 
    const FieldDescr * field_descr,
    int field_index) throw() = 0;

  /// Prepare local array with data to be sent to remote chare for processing
  virtual void prepare_remote (int * n, char ** buffer) throw()
  {};

  /// Accumulate and read data sent from a remote processes
  virtual void update_remote  ( int n, char * buffer) throw()
  {};

  /// Free local array if allocated; NOP if not
  virtual void cleanup_remote (int * n, char ** buffer) throw()
  {};

protected:

  /// Return the filename for the file format and given arguments
  std::string expand_file_name_
  (const std::string * file_name,
   const std::vector<std::string> * file_args) const throw();

private:

  void read_meta_ ( meta_type type, Io * io ) throw();

protected: // attributes

  /// File object for input
  File * file_;

  /// ID of this process
  int process_;

#ifdef CONFIG_USE_CHARM

  /// Counter for ending input
  Counter counter_;

  /// Index of this Input object in Simulation
  size_t index_charm_;

#endif

  /// Simulation cycle for next IO
  int cycle_;

  /// Input counter
  int count_input_;

  /// Simulation time for next IO
  double time_;

  /// Name of the file to read, including format arguments
  std::string file_name_;

  /// Format strings for file name, if any ("cycle", "time", etc.)
  std::vector<std::string> file_args_;

  /// Iterator over field id's
  ItField * it_field_;

  /// I/O Block data accessor
  IoBlock * io_block_;

  /// I/O FieldBlock data accessor
  IoFieldBlock * io_field_block_;

private: // attributes

  /// Only processes with id's divisible by process_stride_ reads
  /// (1: all processes read; 2: 0,2,4,... read; np: root process reads)
  /// Private so that setting must be made through set_process_stride(),

  int process_stride_;
};

#endif /* IO_INPUT_HPP */
