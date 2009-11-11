//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      cello.cpp
 * @brief     Cello main
 * @author    James Bordner
 * @date      Mon Oct  5 15:10:56 PDT 2009
 *
 * $Id$
 *
 *********************************************************************
 */

#include "cello.h"

#include "control.hpp"
#include "monitor.hpp"
#include "parameters.hpp"

void usage(int argc, char ** argv) 
{
  fprintf (stderr,"Usage: %s <filename>\n",argv[0]);
  exit(1);
}

int main(int argc, char ** argv)
{

  try {

    printf ("\n	 ====================================================\n");
    printf ("		      ENZO: THE NEXT GENERATION\n\n");
    printf ("	  A parallel astrophysics and cosmology application\n\n");
    printf ("	     See CELLO_LICENSE for full license agreement\n");
    printf ("	 ====================================================\n\n");

    // INITIALIZE PARALLEL

    // INITALIZE MONITOR

    Monitor monitor;
    monitor.print ("CELLO BEGIN");

    // INITIALIZE CONTROL

    Control control(&monitor);

    // INPUT PARAMETERS

    FILE *fp = 0;
    if (argc == 2) {
      fp = fopen(argv[1],"r");
      if ( !fp ) usage(argc,argv);
    } else {
      usage(argc,argv);
    }

    assert (fp != 0);
  
    control.read_parameters(fp);

    // INITIALIZE SIMULATION

    control.initialize_simulation();

    // RUN SIMULATION

    control.execute_simulation();

    // FINALIZE SIMULATION

    control.terminate_simulation();

    monitor.print ("CELLO END");
  }

  catch (ExceptionBadPointer) {
    printf ("CELLO ERROR: Bad pointer.\n");
  }

}
