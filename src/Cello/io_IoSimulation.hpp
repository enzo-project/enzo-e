// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-24
/// @brief    [\ref Io] Declaration of the IoSimulation class

#ifndef IO_IO_SIMULATION_HPP
#define IO_IO_SIMULATION_HPP

class Simulation;

class IoSimulation : public Io {

  /// @class    IoSimulation
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between Simulation and IO

public: // interface

  /// Constructor
  IoSimulation(const Simulation * simulation) throw();

  /// CHARM++ PUP::able declaration
  PUPable_decl(IoSimulation);

  /// CHARM++ migration constructor
  IoSimulation(CkMigrateMessage *m) : Io(m) {}

  /// Destructor
  virtual ~IoSimulation () throw()
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {

    TRACEPUP;

    // NOTE: change this function whenever attributes change

    Io::pup(p);

    p | rank_; 
    p | cycle_;
    p | time_;
    p | dt_;
  }

  /// Return the ith metadata item associated with the object
  virtual void meta_value 
  (int index, 
   void ** buffer, std::string * name, int * type,
   int * nxd=0, int * nyd=0, int * nzd=0) throw();

  /// Copy the values to the object
  virtual void save_to (void *); 

private: // attributes

  /// Rank of the simulation
  int  rank_; 

  /// Current cycle
  int cycle_;

  /// Current time
  double time_;

  /// Current timestep
  double dt_;

};

#endif /* IO_IO_SIMULATION_HPP */

