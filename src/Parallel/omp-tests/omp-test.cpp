#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "error.hpp"
#include "memory.hpp"
#include "performance.hpp"

#define CODE 1.0 + 3.0*array[k]
#define CODE_FLOPS 2
#define SIZE 128
#define COUNT 500

double gflops (double time)
{
  return CODE_FLOPS*COUNT*SIZE*SIZE*SIZE/time/1024./1024./1024.;
}

int main(int argc, char ** argv)
{
  // Get arguments
  if (argc < 2 || argc > 2) {
    printf ("Usage: %s <num-threads>\n",argv[0]);
    exit(1);
  }
  int threads = atoi(argv[1]);
  // Test argument range
  if (threads < 1 || threads > 32) {
    printf ("threads = %d is out of range\n",threads);
    exit(1);
  }
  // Set Open-MP threads
  omp_set_num_threads(threads);

  double * array = new double [SIZE*SIZE*(SIZE+threads)];

#pragma omp parallel
  {
    int ip_omp = omp_get_thread_num();
    int thread_start        = ip_omp*SIZE/threads;
    int thread_start_padded = ip_omp*SIZE/threads+ip_omp;
    int thread_stop         = (ip_omp+1)*SIZE/threads;
    int thread_stop_padded  = (ip_omp+1)*SIZE/threads+ip_omp;

    Timer t;

    // bring into cache
    for (int iz=0; iz<SIZE; iz++) {
      for (int iy = 0; iy<SIZE; iy++) {
	for (int ix=thread_start; ix<thread_stop; ix++) {
	  int k = ix + SIZE*(iy + iz*SIZE);
	  array[k] = CODE;
	}
      }
    }
    t.start();
    for (int ic = 0; ic < COUNT; ic++) {
      for (int iz=0; iz<SIZE; iz++) {
	for (int iy = 0; iy<SIZE; iy++) {
	  for (int ix=thread_start; ix<thread_stop; ix++) {
	    int k = ix + SIZE*(iy + iz*SIZE);
	    array[k] = CODE;
	  }
	}
      }
    }
    t.stop();
    printf ("X thread-%d (time,gflops) = (%g,%g)  \n",ip_omp,t.value(),
	    gflops(t.value()));

    t.clear();
    for (int iz=thread_start; iz<thread_stop; iz++) {
      for (int iy = 0; iy<SIZE; iy++) {
	for (int ix=0; ix<SIZE; ix++) {
	  int k = ix + SIZE*(iy + iz*SIZE);
	  array[k] = CODE;
	}
      }
    }
    t.start();
    for (int ic = 0; ic < COUNT; ic++) {
      for (int iz=thread_start; iz<thread_stop; iz++) {
	for (int iy = 0; iy<SIZE; iy++) {
	  for (int ix=0; ix<SIZE; ix++) {
	    int k = ix + SIZE*(iy + iz*SIZE);
	    array[k] = CODE;
	  }
	}
      }
    }
    t.stop();
    printf ("Z thread-%d (time,gflops) = (%g,%g)  \n",ip_omp,t.value(),
	    gflops(t.value()));

    t.clear();
    for (int iz=thread_start; iz<thread_stop; iz++) {
      for (int iy = 0; iy<SIZE; iy++) {
	for (int ix=0; ix<SIZE; ix++) {
	  int k = ix + SIZE*(iy + iz*SIZE);
	  array[k] = 1.0 + 3.0*array[k];
	}
      }
    }
    t.start();
    for (int ic = 0; ic < COUNT; ic++) {
      for (int iz=thread_start; iz<thread_stop; iz++) {
	for (int iy = 0; iy<SIZE; iy++) {
	  for (int ix=0; ix<SIZE; ix++) {
	    int k = ix + SIZE*(iy + iz*SIZE);
	    array[k] = 1.0 + 3.0*array[k];
	  }
	}
      }
    }
    t.stop();
    printf ("Z padded thread-%d (time,gflops) = (%g,%g)  \n",ip_omp,t.value(),
	    gflops(t.value()));
  }
}
