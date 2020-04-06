// See LICENSE_CELLO file for license and copyright information

/// @file     test_Face.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Face class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

#include "test_setup_face.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Face");

  // level of block for index 1 and 2

 
  Face * prev_A = nullptr;
  Face * prev_B = nullptr;

  for (int index_test=0; index_test<test::num_face_tests; index_test++) {

    CkPrintf ("Running test %d\n",index_test);
    Index i1 = test::create_index ( test::face_test[index_test].levels_1,
                              test::face_test[index_test].array_1,
                              test::face_test[index_test].tree_1);
    Index i2 = test::create_index ( test::face_test[index_test].levels_2,
                              test::face_test[index_test].array_2,
                              test::face_test[index_test].tree_2);

    const int rank = test::face_test[index_test].rank;
    
    const int ax = test::face_test[index_test].array[0];
    const int ay = test::face_test[index_test].array[1];
    const int az = test::face_test[index_test].array[2];
    const int px = test::face_test[index_test].periodic[0];
    const int py = test::face_test[index_test].periodic[1];
    const int pz = test::face_test[index_test].periodic[2];
    const int ox = test::face_test[index_test].offset[0];
    const int oy = test::face_test[index_test].offset[1];
    const int oz = test::face_test[index_test].offset[2];
    const int NX = test::face_test[index_test].NX;
    const int NY = test::face_test[index_test].NY;
    const int NZ = test::face_test[index_test].NZ;
    const int face = test::face_test[index_test].face;
    const int axis = test::face_test[index_test].axis;
    const int l_degenerate = test::face_test[index_test].l_degenerate;

    const int fx = 0;
    const int fy = 0;
    const int fz = 0;
    
    Face * face_A = new Face (i1,i2,rank,ax,ay,az,px,py,pz,fx,fy,fz);
    Face * face_B = new Face (i2,i1,rank,ax,ay,az,px,py,pz,fx,fy,fz);

    face_A->set_normal(axis,face);
    face_B->set_normal(axis,-face);

    //--------------------------------------------------

    unit_func ("Face()");

    unit_assert (face_A != NULL);
    unit_assert (face_B != NULL);

    unit_assert ((*face_A) == (*face_A));
    unit_assert ((*face_A) == (*face_B)); 
    unit_assert ((*face_B) == (*face_A));
    unit_assert ((*face_B) == (*face_B));
  
    //--------------------------------------------------

    unit_func ("index_block()");

    unit_assert (face_A->index_block() == i1);
    unit_assert (face_B->index_block() == i2);

    if (i1 != i2) {
      unit_assert (face_A->index_block() != i2);
      unit_assert (face_B->index_block() != i1);
    }
  
    //--------------------------------------------------

    unit_func ("index_neighbor()");

    unit_assert (face_A->index_neighbor() == i2);
    unit_assert (face_B->index_neighbor() == i1);
    if (i1 != i2) {
      unit_assert (face_A->index_neighbor() != i1);
      unit_assert (face_B->index_neighbor() != i2);
    }
  
    //--------------------------------------------------

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

    //--------------------------------------------------

    unit_func ("get_adjacency()");

    bool lx,ly,lz;
    face_A->adjacency(&lx,&ly,&lz);

    unit_assert ((lx && axis != 0) || (!lx && axis==0));
    unit_assert ((ly && axis != 1) || (!ly && axis==1));
    unit_assert ((lz && axis != 2) || (!lz && axis==2));
    
    //--------------------------------------------------

    unit_func ("get_offset()");

    int ix0, iy0, iz0;
    face_A->get_offset(&ix0,&iy0,&iz0, NX,NY,NZ);

    if (ix0 != ox*NX/2) {
      CkPrintf ("DEBUG_OFFSET X %d != %d\n",ix0,ox*NX/2);
    }
    if (iy0 != oy*NY/2) {
      CkPrintf ("DEBUG_OFFSET Y %d != %d\n",iy0,oy*NY/2);
    }
    if (iz0 != oz*NZ/2) {
      CkPrintf ("DEBUG_OFFSET Z %d != %d\n",iz0,oz*NZ/2);
    }

    unit_assert (ix0 == ox*NX/2);
    unit_assert (iy0 == oy*NY/2);
    unit_assert (iz0 == oz*NZ/2);
    //--------------------------------------------------

    unit_func ("operator < ()");

    unit_assert ( ((*face_A) < (*face_A)) == false);
    unit_assert ( ! (((*face_A) < (*face_B)) && ((*face_B) < (*face_A))) );
    unit_assert ( ((l_degenerate && ((*face_A) == (*face_B)))
                   || (((*face_A) < (*face_B)) || ((*face_B) < (*face_A)))));
    unit_assert ( ! (((*face_A) < (*face_B)) && ((*face_B) < (*face_A))));

    //--------------------------------------------------

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
    delete prev_A;
    delete prev_B;
    prev_A = face_A;
    prev_B = face_B;
    

  }
  delete prev_A;
  delete prev_B;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

