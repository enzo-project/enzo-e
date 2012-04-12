// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Enzo] Declaration of the EnzoSimulation class

#ifndef ENZO_ENZO_SIMULATION_HPP
#define ENZO_ENZO_SIMULATION_HPP

class EnzoSimulation : public Simulation {

  /// @class    EnzoSimulation
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation base class for Enzo-P

public: // functions

  /// Constructor

#ifdef CONFIG_USE_CHARM

  EnzoSimulation
  ( const char *   parameter_file_name,
    int            n
    ) throw();

#else

  EnzoSimulation
  ( const char *   parameter_file_name,
    const GroupProcess * group_process = 0
    ) throw();
  
#endif


  //==================================================
  // CHARM
  //==================================================

#ifdef CONFIG_USE_CHARM
  /// Initialize an empty EnzoSimulation
  EnzoSimulation();

  /// Initialize a migrated EnzoSimulation
  EnzoSimulation (CkMigrateMessage *m);

  //==================================================

#endif

  /// Destructor
  virtual ~EnzoSimulation() throw();

  /// Override Simulation initialize
  virtual void initialize() throw ();

  /// Finalize the Simulation after running it
  virtual void finalize() throw();

  /// Return an Enzo mesh factory object
  const Factory * factory() const throw();

};

#endif /* ENZO_ENZO_SIMULATION_HPP */

