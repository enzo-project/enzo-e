// See LICENSE_CELLO file for license and copyright information

/// @file     test_Face.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-04-02
/// @brief    Test program for the Face class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

#include "test_setup_face.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  for (int index_test=0; index_test<test::num_face_tests; index_test++) {

    auto test = test::face_test[index_test];

    const int rank = test.rank;
    
    //--------------------------------------------------
    
    unit_func ("Face()");

    int ix = test.face[0];
    int iy = test.face[1];
    int iz = test.face[2];

    const int rx = test.normal[0];
    const int ry = test.normal[1];
    const int rz = test.normal[2];
    
    const int cx = test.child[0];
    const int cy = test.child[1];
    const int cz = test.child[2];

    Face * face_1 = new Face (ix,iy,iz,rx,ry,rz,cx,cy,cz);
    
    if (rx) ix = -ix;
    if (ry) iy = -iy;
    if (rz) iz = -iz;

    Face * face_2 = new Face (ix,iy,iz,rx,ry,rz,cx,cy,cz);

    unit_assert (face_1 != nullptr);
    unit_assert (face_2 != nullptr);
    
    //--------------------------------------------------

    unit_func("face()");

    int ix1,iy1,iz1;
    int ix2,iy2,iz2;

    face_1->get_face(&ix1,&iy1,&iz1);
    face_2->get_face(&ix2,&iy2,&iz2);

    unit_assert(rx ? (ix1 == -ix2) : (ix1 == 0 && ix2 == 0));
    unit_assert(ry ? (iy1 == -iy2) : (iy1 == 0 && iy2 == 0));
    unit_assert(rz ? (iz1 == -iz2) : (iz1 == 0 && iz2 == 0));

    if (rx) unit_assert(std::abs(ix1) + std::abs(ix2) == 2);
    if (ry) unit_assert(std::abs(iy1) + std::abs(iy2) == 2);
    if (rz) unit_assert(std::abs(iz1) + std::abs(iz2) == 2);

    //--------------------------------------------------

    unit_func("normal()");

    int rx1,ry1,rz1;
    int rx2,ry2,rz2;

    face_1->get_normal(&rx1,&ry1,&rz1);
    face_2->get_normal(&rx2,&ry2,&rz2);
    
    unit_assert(rx1 == rx2);
    unit_assert(ry1 == ry2);
    unit_assert(rz1 == rz2);

    //--------------------------------------------------

    unit_func("child()");

    int cx1,cy1,cz1;
    int cx2,cy2,cz2;

    face_1->get_child(&cx1,&cy1,&cz1);
    face_2->get_child(&cx2,&cy2,&cz2);
    
    unit_assert(cx1 == cx2);
    unit_assert(cy1 == cy2);
    unit_assert(cz1 == cz2);

    //--------------------------------------------------

    delete face_1;
    delete face_2;

  }
  
  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

