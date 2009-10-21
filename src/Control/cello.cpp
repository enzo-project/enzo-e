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

#include "monitor.hpp"
#include "parameters.hpp"

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

    // INPUT PARAMETERS

    Parameters parameters;

    if (argc != 2) {
      // Print usage
      char error_message[ERROR_MESSAGE_LENGTH];
      sprintf (error_message,"Usage: %s <filename>\n",argv[0]);
      ERROR_MESSAGE("main",error_message);
    } else {
      // Read in parameters
      FILE * fp = fopen (argv[1],"r");
      if (fp == NULL) {
	char error_message[ERROR_MESSAGE_LENGTH];
	sprintf (error_message,"Filename '%s' unreadable\n",argv[1]);
	ERROR_MESSAGE("main",error_message);
      } else {
	monitor.print ("Reading parameters");
	parameters.read(fp);
      }
    }

    // INITIALIZE SIMULATION

    // RUN SIMULATION

    // FINALIZE SIMULATION

    monitor.print ("CELLO END");
  }

  catch (ExceptionBadPointer) {
    printf ("CELLO ERROR: Bad pointer.\n");
  }

}
