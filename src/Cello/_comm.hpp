// See LICENSE_CELLO file for license and copyright information

/// @file     _comm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-27
/// @brief    Private include file for the \ref Comm component 

#ifndef _COMM_HPP
#define _COMM_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#ifdef CONFIG_USE_CHARM
#  include "pup_stl.h"
#endif

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "mesh_Index.hpp"
#include "comm_CommBlock.hpp"

#endif /* _COMM_HPP */
