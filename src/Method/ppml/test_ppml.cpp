/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
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
/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
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
/** 
 *********************************************************************
 *
 * @file      test_hydro.cpp
 * @brief     Program implementing unit tests for hydrodynamics
 * @author    James Bordner
 * @date      Fri Mar  7 17:11:14 PST 2008
 *
 *********************************************************************
 */

#include "cello.h" 
#include "string.h"
#include "cello_hydro.h"
#include "test_ppml.h"
#include "performance.hpp"

//----------------------------------------------------------------------

void print_usage(const char * name)
{
  printf ("Usage: %s <blast|implosion3> [size] [cycles] [dump-frequency]\n",
	  name);
  exit(1);
}

//----------------------------------------------------------------------

int main(int argc, char * argv[])
{
  enum type_problem problem_type = problem_unknown;

  if (argc < 2) {
    print_usage(argv[0]);
  }

  int argi = 0;

  // Problem type

  if (argc > ++argi) {
    for (int i=0; i<num_problems; i++) {
      if (strcmp(argv[argi],problem_name[i]) == 0) {
	problem_type = type_problem(i);
      }
    }
  }
  if (problem_type == 0) print_usage(argv[0]);

  int size = problem_size[problem_type];
  if (argc > ++argi) size = atoi(argv[argi]);
  int cycles = problem_cycles[problem_type];
  if (argc > ++argi) cycles = atoi(argv[argi]);
  int dump_frequency = 10;
  if (argc > ++argi) dump_frequency = atoi(argv[argi]);

  printf ("problem = %s  size = %d  cycles = %d  dump_frequency = %d\n",
	  problem_name[problem_type], size, cycles, dump_frequency);

  initialize_hydro ();


  switch (problem_type) {
  case problem_blast:
    initialize_ppml(size,cycles);
    break;
  case problem_implosion3:
    initialize_ppml_implosion3(size,cycles);
    break;
  default:
    print_usage(argv[0]);
  }

  float dt;

  int   cycle;
  float time;

  Timer timer;
  timer.start();

  double lower = 0.125*size;
  double upper =   1.0*size;

  for (cycle = 0, time = 0.0;
       (cycle < cycle_stop) && (time < time_stop);
       ++cycle, time += dt) {

    dt =  min (ComputeTimeStep(), time_stop - time);

    SetExternalBoundaryValues();

    if (dump_frequency && (cycle % dump_frequency) == 0) {
      printf ("cycle = %6d seconds = %5.0f sim-time = %6f dt = %6f\n",
	      cycle,timer.value(),time,dt);
      fflush(stdout);
      image_dump(problem_name[problem_type],cycle,lower,upper);
    }

    SolveMHDEquations(cycle, dt);

  }

  printf ("%d %d %g\n",size+6,cycles,timer.value());
  fflush(stdout);

  if (dump_frequency && (cycle % dump_frequency) == 0) {
    SetExternalBoundaryValues();
    image_dump(problem_name[problem_type],cycle,lower,upper);
  }
}

