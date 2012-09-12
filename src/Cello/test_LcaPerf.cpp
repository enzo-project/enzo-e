// See LICENSE_CELLO file for license and copyright information

/// @file     test_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Classname class

#include "main.hpp"
#include "test.hpp"

#include "lcaperf.hpp"

#define N 1000000

float x [N];
float y [N];
float z [N];

//------------------------------------------------------------------------
void foo_flops (LcaPerf * lcaperf,
		int n, float *x, float *y, float *z, const char * region)
//------------------------------------------------------------------------
{
  for (int i=0; i<10; i++) {

    lcaperf->start (region);

    // 20M load
    // 10M flop
    // 10M store
    
    for (int i=0; i<n; i++) {
      x[i] = y[i]*z[i];
    }

    lcaperf->stop (region);

  }
}
#ifdef USE_MPI
//------------------------------------------------------------------------
void foo_mpi (LcaPerf * lcaperf,
	      int count, int size, float *x, float *y, const char * region)
//------------------------------------------------------------------------
{
  int ip,np;
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);

  unit_assert(np > 1);
  unit_assert (0 <= ip && ip < np);

  // Send / Recv / Barrier

  MPI_Status st;
  for (int i=0; i<count; i++) {

    lcaperf->start (region);

    int ip_left  = (ip + np - 1) % np;
    int ip_right = (ip + 1) % np;
    MPI_Send (x, size, MPI_FLOAT, ip_right, i+100*ip,      MPI_COMM_WORLD);
    MPI_Recv (y, size, MPI_FLOAT, ip_left,  i+100*ip_left, MPI_COMM_WORLD,&st);

    MPI_Barrier (MPI_COMM_WORLD);

    lcaperf->stop (region);

  }
}
#endif

//------------------------------------------------------------------------
void foo_enzo(LcaPerf * lcaperf,
	      int n, int count, int size, float *x, float *y, float *z, 
	      const char * region)
//------------------------------------------------------------------------
{

  const int num_levels = 7;
  const int level[] = {0,1,2,2,1,2,2};

  lcaperf->attribute("level",0,LCAP_INT);

  lcaperf->start("EvolveHierarchy");

  for (int i = 0; i < num_levels; i++) {

    lcaperf->attribute("level",&level[i],LCAP_INT);

    lcaperf->start("EvolveLevel");

    foo_flops(lcaperf,n,x,y,z,region);
#ifdef USE_MPI
    foo_mpi(lcaperf,count,size,x,y,region);
#endif

    lcaperf->stop("EvolveLevel");

  }

  lcaperf->attribute("level",0,LCAP_INT);

  lcaperf->stop("EvolveHierarchy");
}

//======================================================================

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("LcaPerf");

  LcaPerf * lcaperf = new LcaPerf;

  
  lcaperf->initialize("out.lcaperf");


  // Initialize arrays

  for (int i=0; i<N; i++) {
    x[i] = 0;
    y[i] = 0;
    z[i] = 0;
  }

  //----------------------------------------------------------------------
  // Built-in counters
  //----------------------------------------------------------------------

  unit_func ("basic");

  lcaperf->begin();

  foo_flops (lcaperf,N,x,y,z,"test_basic");

  lcaperf->end();


  //----------------------------------------------------------------------
  // User counters
  //----------------------------------------------------------------------

  unit_func ("user");

  lcaperf->new_counter("counter_1000",counters_type_absolute);

  lcaperf->begin();

  lcaperf->assign("counter_1000",300);

  foo_flops (lcaperf,N,x,y,z,"test_counter");

  lcaperf->increment("counter_1000",700);

  foo_flops (lcaperf,N,x,y,z,"test_counter");

  unit_assert (lcaperf->value("user", "test_counter", "counter_1000") == 1000);

  lcaperf->end();

  //----------------------------------------------------------------------
  // Attributes
  //----------------------------------------------------------------------

  unit_func ("attributes");

  lcaperf->new_attribute ("level",LCAP_INT);

  lcaperf->begin();

  int level = 1;

  lcaperf->attribute ("level",&level,LCAP_INT);

  foo_flops (lcaperf,N,x,y,z,"test_attribute");

  level = 2;
  lcaperf->attribute ("level",&level,LCAP_INT);

  foo_flops (lcaperf,N,x,y,z,"test_attribute");
  foo_flops (lcaperf,N,x,y,z,"test_attribute");

  int c1 = lcaperf->value("basic", "test_attribute:1", "calls");
  int c2 = lcaperf->value("basic", "test_attribute:2", "calls");

  unit_assert (2*c1 == c2);

  lcaperf->end();

  //----------------------------------------------------------------------
  // PAPI
  //----------------------------------------------------------------------

#ifdef USE_PAPI

  unit_func ("papi");

  lcaperf->begin();

  foo_flops (lcaperf,N,x,y,z,"test_papi");

  double fp = lcaperf->value("papi", "test_papi", "papi-fp-ops");
  //  double tc = lcaperf->value("papi", "test_papi", "papi-tot-cyc");
  // double ld = lcaperf->value("papi", "test_papi", "papi-ld-ins");
  // double sr = lcaperf->value("papi", "test_papi", "papi-sr-ins");

  double eps = 0.01;

  // Ratios: SR=1 LD=2 FP=1
  //  unit_assert (fabs(((sr + fp) - ld)) < (eps * ld));

  unit_assert (fabs(fp - N * 10) < (eps * N*10));
  //  unit_assert (tc > fp);

  //   double tc = lcaperf->value("papi", "test_papi", "papi-tot-cyc");
  // Total cycles > SR + LD + FP
  //  printf ("t%g f%g l%g s%g\n",tc,fp,ld,sr);
  //  unit_assert (tc > fp + ld + sr);
  // NOT ALWAYS TRUE WHEN mpi-off
  // e.g. "t3.90231e+07 f1.0011e+07 l2.00665e+07 s1.00336e+07"

  lcaperf->end();

#endif

  //----------------------------------------------------------------------
  // MPI
  //----------------------------------------------------------------------

#ifdef USE_MPI

  unit_func ("mpi");

  lcaperf->begin();

  foo_mpi (lcaperf,10,100,x,y,"test_mpi");

  int bc = lcaperf->value("mpi", "test_mpi", "mpi-sync-calls");
  int bp = lcaperf->value("mpi", "test_mpi", "mpi-sync-procs");
  int rb = lcaperf->value("mpi", "test_mpi", "mpi-recv-bytes");
  int rc = lcaperf->value("mpi", "test_mpi", "mpi-recv-calls");
  int sb = lcaperf->value("mpi", "test_mpi", "mpi-send-bytes");
  int sc = lcaperf->value("mpi", "test_mpi", "mpi-send-calls");

  unit_assert (bc == 10);
  unit_assert (bp == 10*np);
  unit_assert (rb == 10 * 100 * sizeof(float));
  unit_assert (rc == 10);
  unit_assert (sb == 10 * 100 * sizeof(float));
  unit_assert (sc == 10);

  lcaperf->end();

#endif

  //----------------------------------------------------------------------
  // ENZO print
  //----------------------------------------------------------------------

  unit_func ("mpi");

  lcaperf->new_region ("EvolveHierarchy");
  lcaperf->new_region ("EvolveLevel");
  lcaperf->new_region ("hydro");

  lcaperf->new_attribute ("cycle",LCAP_INT);
  lcaperf->new_attribute ("level",LCAP_INT);

  lcaperf->begin();

  for (int cycle = 100; cycle < 103; cycle ++) {
    
    lcaperf->attribute("cycle",&cycle,LCAP_INT);

    foo_enzo (lcaperf,N,10,100,x,y,z,"hydro");

    lcaperf->print();
  }

  lcaperf->end();

  //----------------------------------------------------------------------
  // attributes key functions
  //----------------------------------------------------------------------

  unit_func("attributes");

  Attributes attributes;
  unit_assert (attributes.keys_match("identical","identical") == true);
  unit_assert (attributes.keys_match("identical","notidentical") == false);

  unit_assert (attributes.keys_match("one:two:three","one:two:three") == true);
  unit_assert (attributes.keys_match("one:three:two","one:two:three") == false);

  unit_assert (attributes.keys_match("one:*:two","one:two:*") == true);
  unit_assert (attributes.keys_match("one:*:*","one:two:*") == true);
  unit_assert (attributes.keys_match("*:*","one:two") == true);
  unit_assert (attributes.keys_match("*:two:four","one:*:three") == false);

  unit_assert (attributes.keys_match("*:*:*","one:two") == false);
  unit_assert (attributes.keys_match("one:two","one:two:three") == false);
  unit_assert (attributes.keys_match("one:two","one:two:*") == false);

  lcaperf->finalize();



  unit_finalize();

  PARALLEL_EXIT;

}

