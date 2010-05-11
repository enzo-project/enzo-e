// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef SIMULATION_SIMULATION_HPP
#define SIMULATION_SIMULATION_HPP

/// @file     simulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-10 16:14:57
/// @brief    Interface file for the Simulation class

#include <vector>
#include "cello.h"
#include "mesh.hpp"
#include "data.hpp"
#include "method.hpp"

class Simulation {

  /// @class    Simulation
  /// @ingroup  Simulation
  /// @brief    Class specifying a simulation to run

public: // interface

  /// Initialize the Simulation object
  Simulation();

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Simulation() throw()
  {
    INCOMPLETE_MESSAGE("Simulation::~Simulation","");
  }

  /// Copy constructor
  Simulation(const Simulation & classname) throw()
  {
    INCOMPLETE_MESSAGE("Simulation::Simulation","");
  }

  /// Assignment operator
  Simulation & operator= (const Simulation & classname) throw()
  {
    INCOMPLETE_MESSAGE("Simulation::operator =","");
    return *this;
  }

  /// Initialize from parameters
  void initialize ();

  /// Return the dimension 1 <= d <= 3, of the simulation
  int dimension() 
  {
    ASSERT ("Simulation::dimension",
	    "lower_ and upper_ vectors are different sizes",
	    lower_.size() == upper_.size()); 
    return lower_.size(); 
  };

  /// Return extents.  Assumes lower[] and upper[] are allocated to at least dimension()
  void domain_extents (int lower[], int upper[])
  {
    if (dimension() >= 1) {
      lower[0] = lower_[0];
      upper[0] = upper_[0];
    }
    if (dimension() >= 2) {
      lower[1] = lower_[1];
      upper[1] = upper_[1];
    }
    if (dimension() >= 3) {
      lower[2] = lower_[2];
      upper[2] = upper_[2];
    }
  }

  /// Return the Mesh object
  Mesh * mesh () { return mesh_; };

  /// Return the Method descriptor
  MethodDescr * methods () { return methods_; };

private: // attributes

  /// Domain extents
  std::vector<double> lower_;
  std::vector<double> upper_;

  /// AMR grid
  Mesh * mesh_;

  /// Methods
  MethodDescr * methods_;

  /// Data descriptions
  DataDescr data_;

};

#endif /* SIMULATION_SIMULATION_HPP */

