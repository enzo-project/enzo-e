/* ------------------------------ DEFINES ------------------------------ */

#ifdef CONFIG_USE_CHARM

#  define PRINTF CkPrintf
#  define ARGC   main->argc
#  define ARGV   main->argv
#  define MAIN  \
      CProxy_Main mainProxy; \
         class Main : public CBase_Main { \
            public:\
               Main(CkArgMsg* main)

#  include "jacobi.decl.h"

#else

#  define PRINTF printf
#  define ARGC   argc
#  define ARGV   argv
#  define MAIN   int main(int argc, char ** argv)

#endif


// ------------------------------ INCLUDES ------------------------------
#ifdef CONFIG_USE_MPI
#   include <mpi.h>
#endif

#include "jacobi_Block.hpp"


