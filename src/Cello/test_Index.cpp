// See LICENSE_CELLO file for license and copyright information

/// @file     test_Index.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Index class

#include "main.hpp"
#include "test.hpp"

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
		     {1,1,0},
		     {1,1,1} };
  index->set_child(1,0,0,0);
  index->set_child(2,1,0,0);
  index->set_child(3,0,1,0);
  index->set_child(4,0,0,1);
  index->set_child(5,0,1,1);
  index->set_child(6,1,1,0);
  index->set_child(7,1,1,1);
  index->set_level(7);
  index->clean();
  unit_func ("set_child()");

  int ic3[3];
  for (int level = 1; level < 7; level++) {

    index->child(level,&ic3[0],&ic3[1],&ic3[2]);
    
    unit_assert(ic3[0] == trace[level-1][0] && 
		ic3[1] == trace[level-1][1] && 
		ic3[2] == trace[level-1][2]);
  }


  int na3[] = { 3,3,3};

  for (int level = 0; level < 7; level++) {

    Index i1,i2;
    i1 = *index;
    i2 = *index;
    int ic3[3];

    i1.set_array(1,0,2);
    i1.set_level(level);
    i1.clean();

    if (level > 0) {

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

    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {

	
	//--------------------------------------------------
 	unit_func("index_neighbor(axis,face)");

	// neighbor is not the same as self
	Index in = i1.index_neighbor(axis,face,na3[axis]);
	unit_assert(in != i1);

	unit_assert(in.level() == level);

	// neighbor's corresponding neighbor is self
	Index in2 = in.index_neighbor(axis,1-face,na3[axis]);
	unit_assert(in2 == i1);

	//--------------------------------------------------
	unit_func("index_neighbor(ix,iy,iz)");

	int i3[3] ;
	i3[0] = 0;
	i3[1] = 0;
	i3[2] = 0;
	i3[axis] = face*2-1;

	// neighbor is not the same as self
	Index ina = i1.index_neighbor(axis,face,na3[axis]);
	Index inb = i1.index_neighbor(i3[0],i3[1],i3[2],na3);
	unit_assert (ina == inb);

	//--------------------------------------------------

	if (level > 0) {

	  unit_func ("index_uncle");
	  Index iu = i1.index_uncle  (axis,face,na3[axis]);

	  unit_assert(i1 != iu);
	  unit_assert(iu.level() == level - 1);

	  unit_func ("index_nibling");

	  i1.child(level,&ic3[0],&ic3[1],&ic3[2]);

	  Index ii = iu.index_nibling(axis,1-face,ic3,na3[axis]);

	  // uncle's corresponding nibling is self
	  unit_assert(i1 == ii); // XXXXX

	}

      }
    }

    unit_func("index_neighbor(ix,iy,iz)");
    for (int ix = -1; ix <= 1; ix++) {
      for (int iy = -1; iy <= 1; iy++) {
	for (int iz = -1; iz <= 1; iz++) {
	  Index in = i1.index_neighbor(ix,iy,iz,na3);
	  if ( ! (ix==0 && iy==0 && iz==0)) unit_assert(in != i1);
	  unit_assert(in.level() == level);
	  Index in2 = in.index_neighbor(-ix,-iy,-iz,na3);
	  unit_assert(in2 == i1);
	}
      }
    }

  }
  

  //	index_uncle (int axis, int face, int narray) const
  // 	index_nibling (int axis, int face, int ic3[3], int narray) const
  //void 	child (int level, int *icx, int *icy, int *icz) const
  // 	child index of this node in parent
  //bool 	is_root () const
  //void 	array (int *ix, int *iy, int *iz) const
  // 	Return the indices of the level-0 node containing this node.
  //int 	level () const
  // 	Return the level of this node.
  //unsigned 	value (int axis) const
  // 	Return the packed bit index for the given axis.
  //void 	set_level (int level)
  // 	Set the level for this node.
  //void 	clean ()
  //void 	set_array (int ix, int iy, int iz)
  // 	Accumulate array part of an index.
  //int 	tree (int axis) const
  // 	Return the packed tree bits for the given axis.
  //void 	set_child (int level, int ix, int iy=0, int iz=0)
  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

