// See LICENSE_CELLO file for license and copyright information

/// @file     test_Adapt.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-10-27
/// @brief    Test program for the Adapt class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"
#include <stdlib.h>

char * black_text (const char c)
{
  char * buffer = new char [80];
  //  sprintf (buffer,"\033[1;30m%c\033[0m",c);
    sprintf (buffer,"%c",c);
  return buffer;
}

char * red_text (const char c)
{
  char * buffer = new char [80];
  sprintf (buffer,"\033[1;31m%c\033[0m",c);
    sprintf (buffer,"%c",c);
  return buffer;
}

char * blue_text (const char c)
{
  char * buffer = new char [80];
  sprintf (buffer,"\033[1;34m%c\033[0m",c);
  sprintf (buffer,"%c",c);
  return buffer;
}

//----------------------------------------------------------------------

bool all_committed(int n, Adapt ** adapt)
{
  for (int i=0; i<n; i++) {
    if (!adapt[i]->is_committed()) return false;
  }
  return true;
}

//----------------------------------------------------------------------

void print_levels (int n, Adapt ** adapt, int min_level, int max_level)
{
  int min=100;
  int max = -100;
  for (int level = max_level+1; level >= min_level-1; level--) {
    for (int i=0; i<n; i++) {
      char c = ' ';
      int level_min = adapt[i]->level_min();
      int level_max = adapt[i]->level_max();
      int level_curr = adapt[i]->level_curr();
      int level_want = adapt[i]->level_want();
      min = std::min(min,level_min);
      max = std::max(max,level_max);
      if (level == level_min && level_min == level_max) {
        c = 'x';
      } else {
        if (level == level_curr) c = '0';
        if (level == level_want) c = '1';
        // if (level_min == level_max) {
        //   if (level == level_min) c = '*';
        // } else {
        //   if (level == level_curr && level_curr == level_want) c = 'o';
        // }
      }
      char * s = nullptr;
      if (level == level_min) s = blue_text(c);
      if (level == level_max) s = red_text(c);
      if (level == level_min && level_min == level_max) s = black_text(c);
      if (s) {
        CkPrintf ((i%2==0) ? "|%s" : "%s",s);
        delete [] s;
      }else
        CkPrintf ((i%2==0) ? "|%c" : "%c",c);
    }
    CkPrintf ("|\n");
  }
  CkPrintf ("C");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->is_committed()?1:0);
  CkPrintf ("\n");
  CkPrintf ("c");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->can_coarsen()?1:0);
  CkPrintf ("\n");
  CkPrintf ("M");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->level_max());
  CkPrintf ("\n");
  CkPrintf ("m");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->level_min());
  CkPrintf ("\n");
  CkPrintf ("c");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->level_curr());
  CkPrintf ("\n");
  CkPrintf ("w");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->level_want());
  CkPrintf ("\n");
  CkPrintf ("n");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->num_neighbors());
  CkPrintf ("\n");
  CkPrintf ("<");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->is_sibling(1)?1:0);
  CkPrintf ("\n");
  CkPrintf (">");
  for (int i0=0; i0<n; i0++) {
    if (adapt[i0]->num_neighbors() >= 3)
      CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->is_sibling(2)?1:0);
    else
      CkPrintf (i0%2?"%c ":"%c",'x');
  }

  CkPrintf ("\n");
}

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Adapt");

  for (int run=0; run<100; run++) {
    
    const int n = 52;
    srand(run);
    std::vector<Index> index;
    index.resize(n);
    ASSERT1("test_Adapt","Index array %d out of range 0 <= i < 1024\n", n,
            (0 <= n && n < 1024));
    for (int i=0; i<n; i++) {
      index[i].set_array(i,0,0);
    }
    int level_curr[n];
    int level_want[n];
    level_curr[0] = 5;
    level_curr[1] = 5;
    for (int i=2; i<n; i+=2) {
      int dx = (rand() % 3) - 1;
      int next = level_curr[i-1] + dx;
      next = std::max(0,next);
      next = std::min(10,next);
      level_curr[i] = next;
      level_curr[i+1] = next;
    }
    for (int i=0; i<n; i++) {
      int dx = (rand() % 3) - 1;
      int want = level_curr[i] + dx;
      want = std::max(0,want);
      want = std::min(10,want);
      level_want[i] = want;
    }
    Adapt * adapt[n];
    for (int i=0; i<n; i++) {
      adapt[i] = new Adapt;
      adapt[i]->set_rank(1);
    }

    // Initialize level bounds for each "block"
    int max_level = -std::numeric_limits<int>::max();
    int min_level = std::numeric_limits<int>::max();
    for (int i0=0; i0<n; i0++) {
      int ip = i0 + 1;
      int im = i0 - 1;
      if (im >= 0 && ip < n) {
        adapt[i0]->allocate_level_bounds(2);
      } else {
        adapt[i0]->allocate_level_bounds(1);
      }
      int k = 0;
      adapt[i0]->set_initial_level_bounds
        (k++,index[i0],level_curr[i0], level_want[i0],true);
      if (im >= 0) {
        adapt[i0]->set_initial_level_bounds
          (k++,index[im],level_curr[im], level_want[im],(i0%2==1));
      }
      if (ip < n) {
        adapt[i0]->set_initial_level_bounds
          (k++,index[ip],level_curr[ip], level_want[ip],(i0%2==0));
      }
      max_level = std::max(max_level,level_curr[i0]);
      max_level = std::max(max_level,level_want[i0]);
      min_level = std::min(min_level,level_curr[i0]);
      min_level = std::min(min_level,level_want[i0]);
    }

    //    print_levels(n,adapt,min_level,max_level);

    unit_func("is_committed()");

    // Update neighbor level bounds

    int iteration = 0;
    while (! all_committed(n,adapt) && iteration < 10) {
      for (int i0=0; i0<n; i0++) {
        int level_min,level_max;
        bool can_coarsen;
        const int im = i0 - 1;
        if (im >= 0) {
          adapt[im]->get_level_bounds(&level_min,&level_max,&can_coarsen);
          adapt[i0]->update_level_bounds(index[im],level_min,level_max,can_coarsen);
        }
        const int ip = i0 + 1;
        if (ip <= n - 1) {
          adapt[ip]->get_level_bounds(&level_min,&level_max,&can_coarsen);
          adapt[i0]->update_level_bounds(index[ip],level_min,level_max,can_coarsen);
        }
        adapt[i0]->evaluate_level_bounds();
      }
      //      CkPrintf ("Iteration %d\n",++iteration);
      //      print_levels(n,adapt,min_level,max_level);
    }

    bool jumps = false;
    bool converged = true;
    bool valid = true;
    bool pairs = true;
    for (int i0=0; i0<n; i0++) {
      int im=std::max(0,i0-1);
      // check for level jumps
      if (im != i0) {
        if (std::abs(adapt[im]->level_min(0) - adapt[i0]->level_min(0)) > 1) jumps=true;
      }
      // check for converged
      if (adapt[i0]->level_min(0) != adapt[i0]->level_max(0)) converged=false;
      // check for refine/coarsen at most once
      if (std::abs(adapt[i0]->level_curr(0) - adapt[i0]->level_min(0)) > 1) valid=false;
      if (std::abs(adapt[i0]->level_curr(0) - adapt[i0]->level_want(0)) > 1) valid=false;
      // check any coarsening done in pairs
      if (adapt[i0]->level_curr(0) > adapt[i0]->level_min(0)) {
        int is = (i0 + 1) - 2*(i0%2);
        if (adapt[is]->level_curr(0) <= adapt[is]->level_min(0)) pairs=false;
      }
      // check coarsening in pairs
    }

    unit_func ("no level jumps");
    unit_assert(jumps == false);
    unit_func ("converged");
    unit_assert(converged == true);
    unit_func ("no temporal jumps");
    unit_assert(valid == true);
    unit_func ("pairs");
    unit_assert(pairs == true);
    if (!pairs || !valid) {
      CkPrintf ("ERROR: C ");
      for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->is_committed()?1:0);
      CkPrintf ("\n");
    }
    for (int i=0; i<n; i++) {
      delete adapt[i];
    }

  }  
  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

