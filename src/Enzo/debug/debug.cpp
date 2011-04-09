// $Id: enzo-p.cpp 2184 2011-04-08 00:58:13Z bordner $
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------

/// @file      enzo-p.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct  5 15:10:56 PDT 2009
/// @todo      support multiple input files
/// @brief     Cello main

//----------------------------------------------------------------------

#define TRACE(MESSAGE) printf ("%s:%d TRACE %s\n",__FILE__,__LINE__,MESSAGE)

#include <stdio.h>

//======================================================================
#include "debug.decl.h"
//======================================================================

#include "mesh_Block.hpp"
#include "enzo_EnzoBlock.hpp"

//----------------------------------------------------------------------


CProxy_Main                proxy_main;

class Main : public CBase_Main {
public:
  Main(CkArgMsg* main)
  {
    CProxy_EnzoBlock test = CProxy_EnzoBlock::ckNew
      (4,4,4,
       0.0, 0.0, 0.0,
       1.0, 1.0, 1.0, 1,
       2,2,2);
    printf ("Blah\n");
  };
};

//======================================================================
#include "debug.def.h"
//======================================================================
