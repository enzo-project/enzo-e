// $Id: parallel_GroupProcess.cpp 1688 2010-08-03 22:34:22Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_GroupProcess.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Aug 11 15:28:49 PDT 2010
/// @brief    Implementation of the GroupProcess class


//----------------------------------------------------------------------

#include "parallel.hpp"

//----------------------------------------------------------------------

GroupProcess * GroupProcess::create (int process_first,
				     int process_last_plus,
				     int process_stride) throw()
{
  GroupProcess * group = 0;

#if defined(CONFIG_USE_MPI)

  group = new GroupProcessMpi
    (process_first, process_last_plus, process_stride);

#elif defined(CONFIG_USE_CHARM)

  INCOMPLETE_MESSAGE ("GroupProcess",
		      "Charm++ not implemented, using serial");

  group = new GroupProcessSerial;

#else

  group = new GroupProcessSerial;

#endif

  ASSERT("GroupProcess","process group creation failed",group != NULL);

  return group;
}
