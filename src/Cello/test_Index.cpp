// See LICENSE_CELLO file for license and copyright information

/// @file     test_Index.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Index class

#include "main.hpp"
#include "test.hpp"
#include "mesh.hpp"

// #include "charm_Index.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Index");

  Index * index = new Index;

  unit_assert (index != NULL);

  int trace[][3] = { {0,0,0},
		     {1,0,0},
		     {0,1,0},
		     {0,0,1},
		     {0,1,1},
		     {1,0,1},
		     {1,1,0},
		     {1,1,1} };

  index->set_child(1,0,0,0);
  index->set_child(2,1,0,0);
  index->set_child(3,0,1,0);
  index->set_child(4,0,0,1);
  index->set_child(5,0,1,1);
  index->set_child(6,1,0,1);
  index->set_child(7,1,1,0);
  index->set_child(8,1,1,1);

  const int max_level = 7;
  index->set_level(max_level);

  index->clean();
  unit_func ("set_child()");

  for (int level = 1; level < max_level; level++) {

    int ic3[3];
    index->child(level,&ic3[0],&ic3[1],&ic3[2]);
    
    unit_assert(ic3[0] == trace[level-1][0] && 
		ic3[1] == trace[level-1][1] && 
		ic3[2] == trace[level-1][2]);
  }


  int na3[] = { 3,3,3};

  for (int level = 0; level < max_level; level++) {

    Index i1,i2;
    i1 = *index;
    i2 = *index;

    i1.set_array(1,0,2);
    i1.set_level(level);
    i1.clean();

    if (level > 0) {

      int ic3[3];
      index->child(level,&ic3[0],&ic3[1],&ic3[2]);
    
      unit_assert(ic3[0] == trace[level-1][0] && 
		  ic3[1] == trace[level-1][1] && 
		  ic3[2] == trace[level-1][2]);


      i2.set_array(1,0,2);
      i2.set_level(level-1);
      i2.clean();

      ic3[0] = trace[level-1][0];
      ic3[1] = trace[level-1][1];
      ic3[2] = trace[level-1][2];

      unit_func("index_parent");
      unit_assert(i1.index_parent() == i2);
      unit_func("index_child");
      unit_assert(i2.index_child(ic3) == i1);

    }

    int f1[3]; // outward face

    for (f1[0]=-1; f1[0]<=1; f1[0]++) {
      for (f1[1]=-1; f1[1]<=1; f1[1]++) {
	for (f1[2]=-1; f1[2]<=1; f1[2]++) {

	  if ( (f1[0]==0 && f1[1]==0 && f1[2]==0) ) continue;

	  //--------------------------------------------------
	  unit_func("index_neighbor()");

	  // neighbor is not the same as self
	  Index n1 = i1.index_neighbor(f1,na3);
	  unit_assert(n1 != i1);

	  unit_assert(n1.level() == level);

	  // neighbor's corresponding neighbor is self
	  int f2[3] = { -f1[0], -f1[1], -f1[2] };
	  Index n2 = n1.index_neighbor(f2, na3);

	  unit_assert(n2 == i1);

	  //--------------------------------------------------

	  if (level > 0) {

	    unit_func ("index_parent");
	    Index p1 = i1.index_parent();
	    unit_assert (p1.level() == i1.level() - 1);

	    unit_func ("index_child");
	    int c1[3]; // child in parent
	    i1.child(level,&c1[0],&c1[1],&c1[2]);
	    Index i2 = p1.index_child(c1);
	    unit_assert (i1 == i2);

	    unit_func ("uncle");

	    int fp[3]; // parent face
	    fp[0] = ((f1[0]==1&&c1[0]==0)||(f1[0]==-1&&c1[0]==1)) ? 0 : f1[0];
	    fp[1] = ((f1[1]==1&&c1[1]==0)||(f1[1]==-1&&c1[1]==1)) ? 0 : f1[1];
	    fp[2] = ((f1[2]==1&&c1[2]==0)||(f1[2]==-1&&c1[2]==1)) ? 0 : f1[2];

	    Index u1 = i1.index_neighbor(f1,na3).index_parent();
	    Index u2 = i1.index_parent().index_neighbor(fp,na3);
	    unit_assert(i1 != u1);
	    unit_assert(u1.level() == level - 1);
	    unit_assert(u1 == u2);

	  }
	}
      }
    }

  }
  


  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

