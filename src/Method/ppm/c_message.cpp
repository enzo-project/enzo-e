/** 
 *********************************************************************
 *
 * @file      c_message.cpp
 * @brief     Display warning and error messages
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @ingroup   Enzo
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */

#include "cello_hydro.h"
 
//----------------------------------------------------------------------
 
int static warning_count = 0;
 
//----------------------------------------------------------------------
 
void c_error (char *sourcefile, Eint32 linenumber)
{
#ifdef USE_MPI
  int  ierr;
#endif
  int error_code;
 
#ifdef USE_MPI
  MPI_Arg id;
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &id);
#else
  int id;
  id = 0;
#endif
 
  printf ("==================\n");
  printf ("=== ENZO ERROR ===   %s: %d   node %d\n",
	  sourcefile,linenumber,id);
  printf ("==================\n");
  fflush(stdout);
 
  error_code = -1;
#ifdef USE_MPI
  ierr = MPI_Abort( MPI_COMM_WORLD, error_code);
#else
  exit(error_code);
#endif
}
 
//----------------------------------------------------------------------
 
void c_warning (char *sourcefile, Eint32 linenumber)
{
#ifdef USE_MPI
  int  ierr;
#endif
 
  ++ warning_count;
 
#ifdef USE_MPI
  MPI_Arg id;
  ierr = MPI_Comm_rank( MPI_COMM_WORLD, &id);
#else
  int id;
  id = 0;
#endif
 
  printf ("--- ENZO WARNING #%d ---   %s: %d   node %d\n",
	  warning_count,sourcefile,linenumber,id);
  fflush(stdout);
 
}

//----------------------------------------------------------------------
extern "C" {
  void FORTRAN_NAME(fc_error) (char *sourcefile, int *linenumber)
  {
    c_error (sourcefile, *linenumber);
  }
}

//----------------------------------------------------------------------

extern "C" {
  void FORTRAN_NAME(fc_warning) (char *sourcefile, int *linenumber)
  {
    c_warning (sourcefile, *linenumber);
  }
}

