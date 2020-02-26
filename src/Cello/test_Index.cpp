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

#define N 8
  Index i8[N+1];

  int trace[][3] = { {0,0,0},
		       {1,0,0},
		       {0,1,0},
		       {0,0,1},
		       {0,1,1},
		       {1,0,1},
		       {1,1,0},
		       {1,1,1} };

  
  // Number of bits (not array size)
  
  int nb3[3] = {3,4,2};
  int na3[3] = {1<<nb3[0], 1<<nb3[1], 1<<nb3[2]};

  for (int i=1; i<N+1; i++) {
    unit_func("operator =()");
    if (i>0) {
      i8[i] = i8[i-1];
      unit_assert(i8[i] == i8[i-1]);
    }
    i8[i].set_level(i);
    i8[i].set_child(i,trace[i-1][0],trace[i-1][1],trace[i-1][2]);
  }

  //==================================================
  // Parent
  //==================================================

  unit_func ("index_parent");

  for (int i=0; i<N; i++) {
    unit_assert (i8[i+1].index_parent() == i8[i]);
  }

  unit_func ("index_parent_level");

  for (int i=0; i<N; i++) {
    for (int l=0; l<=i; l++) {
      unit_assert (i8[i+1].index_ancestor(l) == i8[l]);
    }
  }

  //==================================================
  // Child
  //==================================================

  for (int i=1; i<N; i++) {
    int c3[3];
    i8[i].child(i,c3,c3+1,c3+2);
    unit_func ("child");
    unit_assert (c3[0]==trace[i-1][0] &&
		 c3[1]==trace[i-1][1] &&
		 c3[2]==trace[i-1][2]);
    unit_func ("index_child");
    unit_assert (i8[i].index_child(trace[i][0],trace[i][1],trace[i][2]) == i8[i+1]);
  }

  //==================================================
  // Neighbor
  //==================================================

  unit_func ("index_neighbor");
  for (int axis=0; axis<3; axis++) {
    int w3[3] = {na3[0],na3[1],na3[2]};
    for (int i=0; i<N; i++) {
      int if3[3] = { 0,0,0};

      Index index,index_neighbor;

      // positive direction
      if3[axis] = 1;
      index = i8[i];
      bool l_neighbor = true;
      for (int j=0; j<w3[axis]; j++) {
	index_neighbor = index.index_neighbor(if3,na3);
	l_neighbor = l_neighbor && (index_neighbor != index);
	index = index_neighbor;
      }
      unit_assert (l_neighbor);
      unit_assert (index_neighbor == i8[i]);
      if (index_neighbor != i8[i]) {
	int a3[3];
	int t3[3];
	i8[i].array(a3,a3+1,a3+2);
	i8[i].tree (t3,t3+1,t3+2);
	index_neighbor.array(a3,a3+1,a3+2);
	index_neighbor.tree (t3,t3+1,t3+2);
      }

      // negative direction
      if3[axis] = -1;
      index = i8[i];
      l_neighbor = true;
      for (int j=0; j<w3[axis]; j++) {
	index_neighbor = index.index_neighbor(if3,na3);
	l_neighbor = l_neighbor && (index_neighbor != index);
	index = index_neighbor;
      }
      unit_assert (l_neighbor);
      unit_assert (index_neighbor == i8[i]);
      if (index_neighbor != i8[i]) {
	int a3[3];
	int t3[3];
	i8[i].array(a3,a3+1,a3+2);
	i8[i].tree (t3,t3+1,t3+2);
	index_neighbor.array(a3,a3+1,a3+2);
	index_neighbor.tree (t3,t3+1,t3+2);
      }

      // Update mesh width for next refinement level
      w3[axis] *= 2;
    }
  }


  // ==================================================
  // Array
  // ==================================================

  unit_func ("array");

  int nx = na3[0];
  int ny = na3[1];
  int nz = na3[2];
  Index * index_root = new Index [nx*ny*nz];
  bool l_equal = true;
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i=ix + nx*(iy + ny*iz);
	index_root[i].set_level(0);
	index_root[i].set_array(ix,iy,iz);
	int a3[3];
	index_root[i].array(a3,a3+1,a3+2);
  
	l_equal = l_equal && (a3[0] == ix) && (a3[1] == iy) && (a3[2] == iz);
      }
    }
  }
  unit_assert (l_equal);
  
  // ==================================================
  // Parent (negative levels)
  // ==================================================

  // Level == -1
  int dx=1;
  int dy=nx;
  int dz=nx*ny;
  unit_func ("index_parent (level < 0)");
  l_equal = true;
  for (int iz=0; iz<nz; iz+=2) {
    for (int iy=0; iy<ny; iy+=2) {
      for (int ix=0; ix<nx; ix+=2) {
	int i=ix + nx*(iy + ny*iz);
	Index p = index_root[i].index_parent(-2);
	l_equal = l_equal && (p == index_root[i+dx].index_parent(-2));
	l_equal = l_equal && (p == index_root[i+dy].index_parent(-2));
	l_equal = l_equal && (p == index_root[i+dz].index_parent(-2));
	l_equal = l_equal && (p == index_root[i+dx+dy].index_parent(-2));
	l_equal = l_equal && (p == index_root[i+dy+dz].index_parent(-2));
	l_equal = l_equal && (p == index_root[i+dz+dx].index_parent(-2));
	l_equal = l_equal && (p == index_root[i+dx+dy+dz].index_parent(-2));

	l_equal = l_equal && (p != index_root[i+2*dx].index_parent(-2));
	l_equal = l_equal && (p != index_root[i+2*dy].index_parent(-2));
	l_equal = l_equal && (p != index_root[i+2*dz].index_parent(-2));
	l_equal = l_equal && (p != index_root[i+2*dx+dy].index_parent(-2));
	l_equal = l_equal && (p != index_root[i+2*dx+dz].index_parent(-2));
	l_equal = l_equal && (p != index_root[i+2*dz+dy].index_parent(-2));
      }
    }
  }
  unit_assert (l_equal);

  l_equal = true;
  for (int iz=0; iz<nz; iz+=4) {
    for (int iy=0; iy<ny; iy+=4) {
      for (int ix=0; ix<nx; ix+=4) {
	int i=ix + nx*(iy + ny*iz);
	Index p = index_root[i].index_parent(-2).index_parent(-2);
	l_equal = l_equal && (p == index_root[i+dx].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p == index_root[i+2*dx].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p == index_root[i+3*dx].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p == index_root[i+2*dy].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p == index_root[i+3*dy].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p == index_root[i+2*dz].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p == index_root[i+3*dz].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p != index_root[i+4*dx].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p != index_root[i+4*dy].index_parent(-2).index_parent(-2));
	l_equal = l_equal && (p != index_root[i+4*dz].index_parent(-2).index_parent(-2));
      }
    }
  }
  unit_assert (l_equal);

  //==================================================
  // Neighbor (negative levels)
  //==================================================

  for (int iz=0; iz<nz; iz+=2) {
    for (int iy=0; iy<ny; iy+=2) {
      for (int ix=0; ix<nx; ix+=2) {
	int i=ix + nx*(iy + ny*iz);

	int if3[3] = {1,0,0};
	Index ic = index_root[i];
	Index ip = ic.index_parent(-2);

	// block's neighbor's neighbor's parent
	Index icnnp = ic.index_neighbor(if3,na3);
	icnnp = icnnp.index_neighbor(if3,na3);
	icnnp = icnnp.index_parent(-2);

	// block's parent's neighbor
	Index ipn = ip.index_neighbor(if3,na3);

	unit_assert (icnnp == ipn);
	if (icnnp != ipn) {
	  CkPrintf ("--------------------------------------------------\n");
	  CkPrintf ("   ic    %s\n",ic.bit_string(0,3,nb3).c_str());
	  icnnp = ic.index_neighbor(if3,na3);
	  CkPrintf ("   icn   %s\n",icnnp.bit_string(0,3,nb3).c_str());
	  icnnp = icnnp.index_neighbor(if3,na3);
	  CkPrintf ("   icnn  %s\n",icnnp.bit_string(0,3,nb3).c_str());
	  icnnp = icnnp.index_parent(-2);
	  CkPrintf ("** icnnp %s\n",icnnp.bit_string(0,3,nb3).c_str());
	  CkPrintf ("   ip    %s\n",ip.bit_string(0,3,nb3).c_str());
	  CkPrintf ("** ipn   %s\n",ipn.bit_string(0,3,nb3).c_str());

	}
      }
    }
  }


  //  int na3[] = { 3,3,3};

  Index * index = new Index;
  *index = i8[7];
  for (int level = 0; level < N; level++) {

    Index i1,i2;
    i1 = *index;
    i2 = *index;

    i1.set_array(1,0,2);
    i1.set_level(level);
    //    i1.clean();

    if (level > 0) {

      int ic3[3];
      index->child(level,&ic3[0],&ic3[1],&ic3[2]);
    
      unit_assert(ic3[0] == trace[level-1][0] && 
		  ic3[1] == trace[level-1][1] && 
		  ic3[2] == trace[level-1][2]);


      i2.set_array(1,0,2);
      i2.set_level(level-1);
      //      i2.clean();

      ic3[0] = trace[level-1][0];
      ic3[1] = trace[level-1][1];
      ic3[2] = trace[level-1][2];

      unit_func("index_parent");
      unit_assert(i1.index_parent() == i2);
      unit_func("index_child");
      unit_assert(i2.index_child(ic3) == i1);

    }

    int f1[3]; // outward face

    bool l_neighbor = true;
    bool l_parent = true;
    bool l_child = true;
    bool l_uncle = true;
    
    for (f1[0]=-1; f1[0]<=1; f1[0]++) {
      for (f1[1]=-1; f1[1]<=1; f1[1]++) {
	for (f1[2]=-1; f1[2]<=1; f1[2]++) {

	  if ( (f1[0]==0 && f1[1]==0 && f1[2]==0) ) continue;

	  //--------------------------------------------------
	  unit_func("index_neighbor()");

	  // neighbor is not the same as self
	  Index n1 = i1.index_neighbor(f1,na3);

	  l_neighbor = l_neighbor &&(n1 != i1);
	  l_neighbor = l_neighbor &&(n1.level() == level);

	  // neighbor's corresponding neighbor is self
	  int f2[3] = { -f1[0], -f1[1], -f1[2] };
	  Index n2 = n1.index_neighbor(f2, na3);

	  l_equal = l_equal &&(n2 == i1);

	  //--------------------------------------------------

	  if (level > 0) {

	    unit_func ("index_parent");
	    Index p1 = i1.index_parent();
	    l_parent = l_parent && (p1.level() == i1.level() - 1);

	    unit_func ("index_child");
	    int c1[3]; // child in parent
	    i1.child(level,&c1[0],&c1[1],&c1[2]);
	    Index i2 = p1.index_child(c1);
	    l_child = l_child &&  (i1 == i2);

	    unit_func ("uncle");

	    int fp[3]; // parent face
	    fp[0] = ((f1[0]==1&&c1[0]==0)||(f1[0]==-1&&c1[0]==1)) ? 0 : f1[0];
	    fp[1] = ((f1[1]==1&&c1[1]==0)||(f1[1]==-1&&c1[1]==1)) ? 0 : f1[1];
	    fp[2] = ((f1[2]==1&&c1[2]==0)||(f1[2]==-1&&c1[2]==1)) ? 0 : f1[2];

	    Index u1 = i1.index_neighbor(f1,na3).index_parent();
	    Index u2 = i1.index_parent().index_neighbor(fp,na3);

            l_uncle = l_uncle && (i1 != u1);
            l_uncle = l_uncle && (u1.level() == level - 1);
            l_uncle = l_uncle && (u1 == u2);

	  }
	}
      }
    }

    unit_assert (l_neighbor);
    unit_assert (l_parent);
    unit_assert (l_child);
    unit_assert (l_uncle);
  }

  //==================================================
  // Subtree
  //==================================================

  //  unit_func ("is_in_same_subtree");

  //  unit_assert(false);

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END
