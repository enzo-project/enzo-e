// See LICENSE_CELLO file for license and copyright information

/// @file     io.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-11-01
/// @brief    Include file for the \ref Io component 

#ifndef IO_HPP
#define IO_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include "pngwriter.h"

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "cello.hpp"

#include "error.hpp"
#include "parallel.hpp"
#include "memory.hpp"

#include "simulation.hpp"
#include "mesh.hpp"
#include "disk.hpp"
#include "field.hpp" 

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "io_Io.hpp"
#include "io_IoHierarchy.hpp"
#include "io_IoPatch.hpp"
#include "io_IoBlock.hpp"
#include "io_IoFieldBlock.hpp"

//#include "io_Io.hpp"
// #include "io_IoSimulation.hpp"
// #include "io_IoField.hpp"

#include "io_ItReduce.hpp"
#include "io_ItReduceAvg.hpp"
#include "io_ItReduceSum.hpp"
#include "io_ItReduceMin.hpp"
#include "io_ItReduceMax.hpp"

#include "io_Output.hpp"
#include "io_OutputImage.hpp"
#include "io_OutputData.hpp"

#endif /* IO_HPP */

