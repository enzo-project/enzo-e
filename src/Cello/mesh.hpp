// See LICENSE_CELLO file for license and copyright information

/// @file     mesh.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 11 17:20:03 PST 2010
/// @brief    Include file for the \ref Mesh component 

#ifndef MESH_HPP
#define MESH_HPP

#include "cello.hpp"

#include "charm.hpp"

#include "_performance.hpp"
#include "_monitor.hpp"
#include "_parallel.hpp"
#include "_memory.hpp"
#include "_disk.hpp"
#include "_parameters.hpp"
#include "_control.hpp"
#include "_io.hpp"
#include "_problem.hpp"
#include "_compute.hpp"
#include "_simulation.hpp"
#include "_mesh.hpp"
#include "_data.hpp" 

//----------------------------------------------------------------------
extern void mutex_init_hierarchy();
extern void mutex_init_initial_value();
extern void mutex_init_field_face();
//----------------------------------------------------------------------

#endif /* MESH_HPP */

