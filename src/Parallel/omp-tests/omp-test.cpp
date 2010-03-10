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

//----------------------------------------------------------------------

double gflops (double time)
{
  double size = SIZE*SIZE*SIZE;
  double giga = 1024.0 * 1024.0 * 1024.0;
  double flops = COUNT*CODE_FLOPS*size;

  return (flops/giga)/time;
}

//----------------------------------------------------------------------

void update_y(Timer & timer,
	      int thread_start,
	      int thread_stop,
	      double * array,
	      double * sum)
{
    *sum = 0;
    for (int iz = 0; iz<SIZE; iz++) {
      for (int iy=thread_start; iy<thread_stop; iy++) {
	for (int ix=0; ix<SIZE; ix++) {
	  int k = ix + SIZE*(iy + iz*SIZE);
	  array[k] = k;
	  *sum += array[k];
	}
      }
    }
    timer.start();
    for (int ic = 0; ic < COUNT; ic++) {
      for (int iz = 0; iz<SIZE; iz++) {
	for (int iy=thread_start; iy<thread_stop; iy++) {
	  for (int ix=0; ix<SIZE; ix++) {
	    int k = ix + SIZE*(iy + iz*SIZE);
	    array[k] = CODE;
	  }
	}
      }
    }
    timer.stop();
}


void print_time (Timer & timer,char axis,int ip_omp,double sum)
{
  printf ("(axis,thread,time,gflops,sum) = (%c,%d,%g,%g,%g)\n",
	  axis,ip_omp,timer.value(), gflops(timer.value()),sum);
}

//----------------------------------------------------------------------

int main(int argc, char ** argv)
{
  // Parse arguments

  if (argc < 2 || argc > 2) {
    printf ("Usage: %s <num-threads>\n",argv[0]);
    exit(1);
  }

  int threads = atoi(argv[1]);

  // Verify argument range

  if (threads < 1 || threads > 32) {
    printf ("threads = %d is out of range\n",threads);
    exit(1);
  }

  // Set Open-MP threads

  omp_set_num_threads(threads);

  double * array = new double [SIZE*SIZE*(SIZE+threads)];

#pragma omp parallel
  {
    int ip_omp              = omp_get_thread_num();
    int thread_start        = ip_omp*SIZE/threads;
    int thread_stop         = (ip_omp+1)*SIZE/threads;
    double sum = 0;

    Timer timer;

    // bring into cache
    sum = 0;
    for (int iz=0; iz<SIZE; iz++) {
      for (int iy = 0; iy<SIZE; iy++) {
	for (int ix=thread_start; ix<thread_stop; ix++) {
	  int k = ix + SIZE*(iy + iz*SIZE);
	  array[k] = k;
	  sum += array[k];
	}
      }
    }
    timer.start();
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
    timer.stop();
    printf ("(axis,thread,time,gflops,sum) = (X,%d,%g,%g,%g)\n",
	    ip_omp,timer.value(), gflops(timer.value()),sum);


    update_y(timer,thread_start,thread_stop,array,&sum);
    print_time (timer,'Y',ip_omp,sum);

    timer.clear();
    sum = 0;
    for (int iz=thread_start; iz<thread_stop; iz++) {
      for (int iy = 0; iy<SIZE; iy++) {
	for (int ix=0; ix<SIZE; ix++) {
	  int k = ix + SIZE*(iy + iz*SIZE);
	  array[k] = k;
	  sum += array[k];
	}
      }
    }
    timer.start();
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
    timer.stop();
    printf ("(axis,thread,time,gflops,sum) = (Z,%d,%g,%g,%g)\n",
	    ip_omp,timer.value(), gflops(timer.value()),sum);
  }
}
