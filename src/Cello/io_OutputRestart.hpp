// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputRestart.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2011-09-26
/// @brief    [\ref Io] Declaration for the OutputRestart class

#ifndef IO_OUTPUT_RESTART_HPP
#define IO_OUTPUT_RESTART_HPP

class ItField;
class Parameters;

class OutputRestart : public OutputData {

  /// @class    OutputRestart
  /// @ingroup  Io
  /// @brief [\ref Io] class for reading and writing Restart files

public: // functions

  /// Create an uninitialized OutputRestart object
  OutputRestart(const Factory * factory, Parameters * parameters) throw();

  ///
  ~OutputRestart() throw()
  {}

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change

    OutputData::pup(p);

    p | param_name_;
    p | param_args_;
  }

#endif

public: // virtual functions

  /// Finalize output
  virtual void finalize () throw ();

  /// Write Simulation data to disk
  virtual void write ( const Simulation * simulation ) throw();

private:

  /// Name of the parameter file to write, including format arguments 
  std::string param_name_;

  /// Format strings for param_name_, if any ("cycle", "time", etc.)
  std::vector<std::string> param_args_;

};

#endif /* IO_OUTPUT_RESTART_HPP */
