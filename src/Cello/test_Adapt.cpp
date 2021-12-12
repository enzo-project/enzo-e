#include "main.hpp"
#include "test.hpp"
#include "mesh.hpp"
#include <stdlib.h>

#define NEW_ADAPT_TESTS
//----------------------------------------------------------------------

bool all_converged(int n, Adapt ** adapt)
{
  for (int i=0; i<n; i++) {
    if (!adapt[i]->is_converged()) return false;
  }
  return true;
}

//----------------------------------------------------------------------

void print_levels (int n, Adapt ** adapt, int min_level, int max_level,
                   int * level_curr, int * level_want)
{
  int min=100;
  int max = -100;
  for (int level = max_level+1; level >= min_level-1; level--) {
    for (int i=0; i<n; i++) {
      char c = ' ';
      int level_min = adapt[i]->level_min();
      int level_max = adapt[i]->level_max();
      min = std::min(min,level_min);
      max = std::max(max,level_max);
      if (level == level_min && level_min == level_max) {
        c = 'x';
      } else {
        if (level == level_curr[i]) c = '0';
        if (level == level_want[i]) c = '1';
        // if (level_min == level_max) {
        //   if (level == level_min) c = '*';
        // } else {
        //   if (level == level_curr && level_curr == level_want) c = 'o';
        // }
      }
      CkPrintf ((i%2==0) ? "|%c" : "%c",c);
    }
    CkPrintf ("|\n");
  }
  CkPrintf ("C");
  for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->is_converged()?1:0);
  CkPrintf ("\n");
}

//----------------------------------------------------------------------

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Adapt");

  const int num_runs = 100;
  const int run_length = 52;
  const int MAX_LEVEL = 100;

  for (int run=0; run<num_runs; run++) {

    const int n = run_length;
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
    index[0].set_level(5);
    index[1].set_level(5);
    for (int i=2; i<n; i+=2) {
      int dx = (rand() % 3) - 1;
      const int curr = level_curr[i-1];
      int next = curr + dx;
      next = std::max(0,next);
      next = std::min(10,next);
      level_curr[i] = next;
      level_curr[i+1] = next;
      index[i].set_level(next);
      index[i+1].set_level(next);
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
      adapt[i]->set_max_level(MAX_LEVEL);
      adapt[i]->set_index(index[i]);
    }

    // Initialize level bounds for each "block"
    int max_level = -std::numeric_limits<int>::max();
    int min_level = std::numeric_limits<int>::max();
    for (int i0=0; i0<n; i0++) {
      int ip = i0 + 1;
      int im = i0 - 1;
      // initialize indices
      // initialize level bounds
      int min = level_want[i0];
      int now = level_curr[i0];
      int max = level_curr[i0] + 1;
      adapt[i0]->initialize_self (index[i0],min, now);
      if (im >= 0) {
        adapt[i0]->insert_neighbor (index[im],(i0%2==1));
      }

      if (ip < n) {
        adapt[i0]->insert_neighbor (index[ip],(i0%2==0));
      }
      max_level = std::max(max_level,level_curr[i0]);
      max_level = std::max(max_level,level_want[i0]);
      min_level = std::min(min_level,level_curr[i0]);
      min_level = std::min(min_level,level_want[i0]);
    }

    print_levels(n,adapt,min_level,max_level,level_curr,level_want);

    unit_func("is_converged()");

    // Update neighbor level bounds

    int iteration = 0;
    while (! all_converged(n,adapt) && iteration < 10) {
      for (int i0=0; i0<n; i0++) {
        int level_min,level_max;
        bool can_coarsen;
        const int im = i0 - 1;
        if (im >= 0) {
          adapt[im]->get_level_bounds(&level_min,&level_max,&can_coarsen);
          adapt[i0]->update_neighbor(index[im],level_min,level_max,can_coarsen);
        }
        const int ip = i0 + 1;
        if (ip <= n - 1) {
          adapt[ip]->get_level_bounds(&level_min,&level_max,&can_coarsen);
          adapt[i0]->update_neighbor(index[ip],level_min,level_max,can_coarsen);
        }
        adapt[i0]->update_bounds();
      }
      CkPrintf ("Iteration %d\n",++iteration);
      fflush(stdout);
      print_levels(n,adapt,min_level,max_level,level_curr,level_want);
    }

    bool jumps = false;
    bool converged = true;
    bool valid = true;
    bool pairs = true;
    for (int i0=0; i0<n; i0++) {
      int im=std::max(0,i0-1);
      // check for level jumps
      if (im != i0) {
        if (std::abs(adapt[im]->level_min() - adapt[i0]->level_min()) > 1) jumps=true;
      }
      // check for converged
      if (adapt[i0]->level_min() != adapt[i0]->level_max()) converged=false;
      // check for refine/coarsen at most once
      if (std::abs(level_curr[i0] - adapt[i0]->level_min()) > 1) valid=false;
      if (std::abs(level_curr[i0] - level_want[i0]) > 1) valid=false;
      // check any coarsening done in pairs
      if (level_curr[i0] > adapt[i0]->level_min()) {
        int is = (i0 + 1) - 2*(i0%2);
        if (level_curr[is] <= adapt[is]->level_min()) pairs=false;
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
      for (int i0=0; i0<n; i0++) CkPrintf (i0%2?"%1X ":"%1X",adapt[i0]->is_converged()?1:0);
      CkPrintf ("\n");
    }
    for (int i=0; i<n; i++) {
      delete adapt[i];
    }

  }

  // Test refine / coarsen neighbors

  //   +-------+-------+-------+
  //   |       |       |       |
  //   |   6   |   7   |   8   |       +---+---+---+---+---+
  //   |       |       |       |       |   |   |   |   |   |
  //   +-------+-------+-------+       +---+---+---+---+---+
  //   |       |       |       |  >>>  |   |   |       |   |
  //   |   3   |   4   |   5   |       +---+---+   4   +---+
  //   |       |       |       |  <<<  |   |   |       |   |
  //   +-------+-------+-------+       +---+---+---+---+---+
  //   |       |       |       |       |   |   |   |   |   |
  //   |   0   |   1   |   2   |       +---+---+---+---+---+
  //   |       |       |       |       |   |   |   |   |   |
  //   +-------+-------+-------+       +---+---+---+---+---+


  Index index_refine_start[9];

  index_refine_start[0].set_array(0,0,0);
  index_refine_start[0].push_child(0,0,0);
  index_refine_start[1].set_array(0,0,0);
  index_refine_start[1].push_child(1,0,0);
  index_refine_start[2].set_array(1,0,0);
  index_refine_start[2].push_child(0,0,0);
  index_refine_start[3].set_array(0,0,0);
  index_refine_start[3].push_child(0,1,0);
  index_refine_start[4].set_array(0,0,0);
  index_refine_start[4].push_child(1,1,0);
  index_refine_start[5].set_array(1,0,0);
  index_refine_start[5].push_child(0,1,0);
  index_refine_start[6].set_array(0,1,0);
  index_refine_start[6].push_child(0,0,0);
  index_refine_start[7].set_array(0,1,0);
  index_refine_start[7].push_child(1,0,0);
  index_refine_start[8].set_array(1,1,0);
  index_refine_start[8].push_child(0,0,0);

  Adapt adapt;

  adapt.set_rank(2);
  adapt.set_index(index_refine_start[4]);
  adapt.insert_neighbor(index_refine_start[0]);
  adapt.insert_neighbor(index_refine_start[1]);
  adapt.insert_neighbor(index_refine_start[2]);
  adapt.insert_neighbor(index_refine_start[3]);
  adapt.insert_neighbor(index_refine_start[5]);
  adapt.insert_neighbor(index_refine_start[6]);
  adapt.insert_neighbor(index_refine_start[7]);
  adapt.insert_neighbor(index_refine_start[8]);

  adapt.print("ADAPT_START");

  adapt.refine_neighbor(index_refine_start[0]);
  adapt.refine_neighbor(index_refine_start[1]);
  adapt.refine_neighbor(index_refine_start[2]);
  adapt.refine_neighbor(index_refine_start[3]);
  adapt.refine_neighbor(index_refine_start[5]);
  adapt.refine_neighbor(index_refine_start[6]);
  adapt.refine_neighbor(index_refine_start[7]);
  adapt.refine_neighbor(index_refine_start[8]);
  
  // Test refine / coarsen self

  //   +-------+-------+-------+        +-------+-------+-------+
  //   |       |       |       |        |       |       |       |
  //   |   6   |   7   |   8   |        |   6   |   7   |   8   |
  //   |       |       |       |        |       |       |       |
  //   +-------+-------+-------+        +-------+-------+-------+
  //   |       |       |       |  >>>   |       | B | C |       |
  //   |   3   |   4   |   5   |        |   3   +---+---+   5   |
  //   |       |       |       |  <<<   |       | 9 | A |       |
  //   +-------+-------+-------+        +-------+-------+-------+
  //   |       |       |       |        |       |       |       |
  //   |   0   |   1   |   2   |        |   0   |   1   |   2   |
  //   |       |       |       |        |       |       |       |
  //   +-------+-------+-------+        +-------+-------+-------+

  adapt.print("ADAPT_REFINE");

  const int n = adapt.num_neighbors();
  Index * index_list = new Index[n];
  for (int i=0; i<n; i++) index_list[i] = adapt.index(i);
  for (int i=0; i<n; i++) {
    adapt.coarsen_neighbor(index_list[i]);
  }
  adapt.print("ADAPT");

  Adapt adapt_refine[4];

  int ic3[3] = {0,0,0};
  ic3[0] = 0; ic3[1] = 0;
  adapt_refine[0].set_rank(2);
  adapt_refine[0].refine(adapt, ic3);
  adapt_refine[0].print("REFINE-00");
  ic3[0] = 1; ic3[1] = 0;
  adapt_refine[1].set_rank(2);
  adapt_refine[1].refine(adapt, ic3);
  adapt_refine[1].print("REFINE-10");
  ic3[0] = 0; ic3[1] = 1;
  adapt_refine[2].set_rank(2);
  adapt_refine[2].refine(adapt, ic3);
  adapt_refine[2].print("REFINE-01");
  ic3[0] = 1; ic3[1] = 1;
  adapt_refine[3].set_rank(2);
  adapt_refine[3].refine(adapt, ic3);
  adapt_refine[3].print("REFINE-11");

  Adapt adapt_coarsen;

  adapt_coarsen.set_rank(2);

  adapt_coarsen.coarsen(adapt_refine[0]);
  adapt_coarsen.print("COARSEN-00");
  adapt_coarsen.coarsen(adapt_refine[1]);
  adapt_coarsen.print("COARSEN-10");
  adapt_coarsen.coarsen(adapt_refine[2]);
  adapt_coarsen.print("COARSEN-01");
  adapt_coarsen.coarsen(adapt_refine[3]);
  adapt_coarsen.print("COARSEN-11");

  adapt_coarsen.print("COARSEN");
  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

