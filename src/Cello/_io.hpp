// See LICENSE_CELLO file for license and copyright information

/// @file     _io.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-11-01
/// @brief    Private include file for the \ref Io component 

#ifndef _IO_HPP
#define _IO_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <limits>
#include "pngwriter.h"

//----------------------------------------------------------------------
// Typedefs
//----------------------------------------------------------------------

enum meta_type {
  meta_type_file,
  meta_type_group
};

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "io_Io.hpp"
#include "io_IoHierarchy.hpp"
#include "io_IoPatch.hpp"
#include "io_IoBlock.hpp"
#include "io_IoFieldBlock.hpp"

#include "io_ItReduce.hpp"
#include "io_ItReduceAvg.hpp"
#include "io_ItReduceSum.hpp"
#include "io_ItReduceMin.hpp"
#include "io_ItReduceMax.hpp"

#include "io_Input.hpp"
#include "io_InputData.hpp"

#include "io_Output.hpp"
#include "io_OutputImage.hpp"
#include "io_OutputData.hpp"
#include "io_OutputRestart.hpp"

#include "io_Schedule.hpp"


#endif /* _IO_HPP */

