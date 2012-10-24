// See LICENSE_CELLO file for license and copyright information

/// @file     _performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @brief    Private include file for the \ref Performance component

#ifndef _PERFORMANCE_HPP
#define _PERFORMANCE_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <vector>
#include <map>
#include <stack>
#include <string>
#include <sstream>
#include <sys/resource.h>

#ifdef __linux__
#   include <unistd.h>
#endif
#ifdef CONFIG_USE_PAPI
#  include "papi.h"
#endif

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "performance_Timer.hpp"
#include "performance_PerfCounters.hpp"
#include "performance_Papi.hpp"
#include "performance_Performance.hpp"


#endif /* _PERFORMANCE_HPP */
