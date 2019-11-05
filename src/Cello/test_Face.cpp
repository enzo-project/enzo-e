// See LICENSE_CELLO file for license and copyright information

/// @file     test_Face.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Face class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

Index create_index (const int levels, int array[3], int tree[][3])
{
  Index index;
  index.clear();
  index.set_array (array[0],array[1],array[2]);
  for (int i=0; i<levels; i++) {
    index.set_child(i+1,tree[i][0],tree[i][1],tree[i][2]);
  }
  index.set_level(levels);
  return index;
};

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Face");

  // level of block for index 1 and 2

#include "test_setup_face.hpp"
  
  for (int index_test=0; index_test<num_tests; index_test++) {

    CkPrintf ("Running test %d\n",index_test);
    Index i1 = create_index ( test[index_test].levels_1,
                              test[index_test].array_1,
                              test[index_test].tree_1);
    Index i2 = create_index ( test[index_test].levels_2,
                              test[index_test].array_2,
                              test[index_test].tree_2);

    const int rank = test[index_test].rank;
    
    const int ax = test[index_test].array[0];
    const int ay = test[index_test].array[1];
    const int az = test[index_test].array[2];
    const int px = test[index_test].periodic[0];
    const int py = test[index_test].periodic[1];
    const int pz = test[index_test].periodic[2];
    const int NX = test[index_test].NX;
    const int NY = test[index_test].NY;
    const int NZ = test[index_test].NZ;
    const int face = test[index_test].face;
    const int axis = test[index_test].axis;
    const int l_degenerate = test[index_test].l_degenerate;

    const int fx = 0;
    const int fy = 0;
    const int fz = 0;
    
    Face * face_A = new Face (i1,i2,rank,ax,ay,az,px,py,pz,fx,fy,fz);
    Face * face_B = new Face (i2,i1,rank,ax,ay,az,px,py,pz,fx,fy,fz);

    face_A->set_normal(axis,face);
    face_B->set_normal(axis,-face);

    unit_func ("Face()");

    unit_assert (face_A != NULL);
    unit_assert (face_B != NULL);

    unit_assert ((*face_A) == (*face_A));
    unit_assert ((*face_A) == (*face_B)); 
    unit_assert ((*face_B) == (*face_A));
    unit_assert ((*face_B) == (*face_B));
  
    unit_func ("index_block()");

    unit_assert (face_A->index_block() == i1);
    unit_assert (face_B->index_block() == i2);

    if (i1 != i2) {
      unit_assert (face_A->index_block() != i2);
      unit_assert (face_B->index_block() != i1);
    }
  
    unit_func ("index_neighbor()");

    unit_assert (face_A->index_neighbor() == i2);
    unit_assert (face_B->index_neighbor() == i1);
    if (i1 != i2) {
      unit_assert (face_A->index_neighbor() != i1);
      unit_assert (face_B->index_neighbor() != i2);
    }
  
    unit_func ("get_subface()");

    int f1x,f1y,f1z;
    int f2x,f2y,f2z;

    face_A->get_subface(&f1x,&f1y,&f1z);
    face_B->get_subface(&f2x,&f2y,&f2z);

    unit_assert (f1x == 0);
    unit_assert (f1y == 0);
    unit_assert (f1z == 0);
    unit_assert (f2x == 0);
    unit_assert (f2y == 0);
    unit_assert (f2z == 0);

    unit_func ("get_dimensions()");

    int nx,ny,nz;
    face_A->get_dimensions(&nx,&ny,&nz,
                          NX,NY,NZ);
    const float * size3 = test[index_test].size;

    int nx_exp=(NX&&size3[0]==0) ? 1 : ( (size3[0]!=-1) ? size3[0]*NX : 0);
    int ny_exp=(NY&&size3[1]==0) ? 1 : ( (size3[1]!=-1) ? size3[1]*NY : 0);
    int nz_exp=(NZ&&size3[2]==0) ? 1 : ( (size3[2]!=-1) ? size3[2]*NZ : 0);

    CkPrintf ("DEBUG Testing %d == %d  %d == %d  %d == %d\n",
              nx,nx_exp,ny, ny_exp, nz, nz_exp);
    
    unit_assert (nx == nx_exp);
    unit_assert (ny == ny_exp);
    unit_assert (nz == nz_exp);
    
    unit_func ("operator < ()");

    unit_assert ( ((*face_A) < (*face_A)) == false);
    unit_assert ( ! (((*face_A) < (*face_B)) && ((*face_B) < (*face_A))) );
    unit_assert ( (l_degenerate && ((*face_A) == (*face_B))
                   || (((*face_A) < (*face_B)) || ((*face_B) < (*face_A)))));
    unit_assert ( ! (((*face_A) < (*face_B)) && ((*face_B) < (*face_A))));

    unit_func ("operator == ()");

    unit_assert ( ((*face_A) == (*face_A)));
    unit_assert ( ((*face_A) == (*face_B)));
    unit_assert ( ((*face_B) == (*face_A)));
    if (prev_A != nullptr) {
      unit_assert ( ! ((*prev_A) == (*face_A)));
      unit_assert ( ! ((*prev_A) == (*face_B)));
      unit_assert ( ! ((*face_A) == (*prev_A)));
      unit_assert ( ! ((*face_B) == (*prev_A)));
    }
    if (prev_B != nullptr) {
      unit_assert ( ! ((*prev_B) == (*face_A)));
      unit_assert ( ! ((*prev_B) == (*face_B)));
      unit_assert ( ! ((*face_A) == (*prev_B)));
      unit_assert ( ! ((*face_B) == (*prev_B)));
    }

    delete face_A;
    delete face_B;
  }

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

