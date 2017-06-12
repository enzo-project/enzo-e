// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputCheckpoint.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2011-09-26
/// @brief    [\ref Io] Declaration for the OutputCheckpoint class

#ifndef IO_OUTPUT_CHECKPOINT_HPP
#define IO_OUTPUT_CHECKPOINT_HPP

class OutputCheckpoint : public Output {

  /// @class    OutputCheckpoint
  /// @ingroup  Io
  /// @brief [\ref Io] class for reading and writing Checkpoint files

public: // functions

  /// Empty constructor for Charm++ pup()
  OutputCheckpoint() throw() { }

  /// Create an uninitialized OutputCheckpoint object
  OutputCheckpoint(int index, 
		   const Factory * factory, 
		   const FieldDescr * field_descr,
		   const ParticleDescr * particle_descr,
		   Config * config, 
		   int process_count) throw();

  /// Destructor
  ~OutputCheckpoint() throw()  { }

  /// Charm++ PUP::able declarations
  PUPable_decl(OutputCheckpoint);

  /// Charm++ PUP::able migration constructor
  OutputCheckpoint (CkMigrateMessage *m) : Output (m) { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

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
  virtual void write_field_data
  ( const FieldData * field_data, 
    const FieldDescr * field_descr,
    int index_field) throw()
  { /* EMPTY */ }

  /// Write local particle to disk
  virtual void write_particle_data
  ( const ParticleData * particle_data, 
    const ParticleDescr * particle_descr,
    int index_particle) throw()
  { /* EMPTY */ }

private: // private functions

  /// Read the restart_file_ and update Simulation::config() with
  /// updated values
  void update_config_();

  private: // attributes

  /// Name of parameter file to read on restart for updated parameters
  std::string restart_file_;

};

#endif /* IO_OUTPUT_CHECKPOINT_HPP */
