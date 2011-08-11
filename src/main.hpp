// See LICENSE_CELLO file for license and copyright information

/// @file     main.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file main.cpp
///
/// Detailed description of file main.cpp

#ifdef CONFIG_USE_CHARM

#include "cello.hpp"

#include "parallel.hpp"
#include "monitor.hpp"

class Factory;
#include "main.decl.h"
extern CProxy_Main proxy_main;

//----------------------------------------------------------------------

class Main : public CBase_Main {

public:

  /// Initialize the Main chare (defined in the calling program)
  Main(CkArgMsg* main);
  
  /// Exit the program
  void p_exit(int count);

private:

   int count_exit_; 
   Monitor * monitor_;

};
#endif
