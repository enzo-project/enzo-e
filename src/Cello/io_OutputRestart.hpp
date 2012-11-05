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

  /// Empty constructor for Charm++ pup()
  OutputRestart() throw() { }

  /// Create an uninitialized OutputRestart object
  OutputRestart(int index, 
		const Factory * factory, 
		Config * config, 
		int process_count) throw();

  /// Destructor
  ~OutputRestart() throw()  { }

#ifdef CONFIG_USE_CHARM

  /// Charm++ PUP::able declarations
  PUPable_decl(OutputRestart);

  /// Charm++ PUP::able migration constructor
  OutputRestart (CkMigrateMessage *m) : Output (m) { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#endif

public: // virtual functions

  /// Open (or create) a file for IO
  virtual void open () throw()
  { /* EMPTY */ };

  /// Close file for IO
  virtual void close () throw()
  { /* EMPTY */ };
  
  /// Write Simulation data to disk
  virtual void write_simulation ( const Simulation * simulation ) throw();

  /// Write local field to disk
  virtual void write_field_block
  ( const FieldBlock * field_block, 
    const FieldDescr * field_descr,
    int field_index) throw()
  { /* EMPTY */ }

private:

  /// Name of the parameter file to write, including format arguments 
  std::string dir_name_;

  /// Format strings for dir_name_, if any ("cycle", "time", etc.)
  std::vector<std::string> dir_args_;

};

#endif /* IO_OUTPUT_RESTART_HPP */
