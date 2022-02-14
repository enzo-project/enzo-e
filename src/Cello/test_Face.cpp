// See LICENSE_CELLO file for license and copyright information

/// @file     test_Face.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2020-04-02
/// @brief    Test program for the Face class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("FaceFluxes");

  struct face_test_type {
    int child[3];
    int size[3];
    int centered[3];
    double cell_width[3];
    double time_step_1;
    double time_step_2;
  } test =
      {
       { 0, 0, 0},
       {10, 6, 8 },
       { 0, 0, 0},
       { 0.25, 0.5, 1.0},
       0.125, 0.125
      };

  const int num_level = 3;
  const int level_1[] = {1,2,3};
  const int level_2[] = {1,3,2};
  
  for (int rank = 2; rank <=3; ++rank) {
    for (int iaxis=0; iaxis<rank; iaxis++) {
      for (int iface=0; iface<2; ++iface) {
        for (int ilevel=0; ilevel<num_level; ++ilevel) {
        
          const int L_1 = level_1[ilevel];
          const int L_2 = level_2[ilevel];
          double h3_1[3];
          h3_1[0] = test.cell_width[0];
          h3_1[1] = (rank >= 2) ? test.cell_width[1] : 0;
          h3_1[2] = (rank >= 3) ? test.cell_width[2] : 0;
          double h3_2[3] = {h3_1[0],h3_1[1],h3_1[2]};
          
          if (L_1 < L_2) { h3_2[0]*=0.5; h3_2[1]*=0.5; h3_2[2]*=0.5;}
          if (L_1 > L_2) { h3_2[0]*=2.0; h3_2[1]*=2.0; h3_2[2]*=2.0;}
          int fx = iaxis==0 ? iface*2-1 : 0;
          int fy = iaxis==1 ? iface*2-1 : 0;
          int fz = iaxis==2 ? iface*2-1 : 0;
          
          Face * face_1 = new Face( fx, fy, fz,iaxis,iface);
          Face * face_2 = new Face(-fx,-fy,-fz,iaxis,iface);
          //--------------------------------------------------
    
          unit_func ("Face()");
    
          unit_assert (face_1 != nullptr);
          unit_assert (face_2 != nullptr);
    
          //--------------------------------------------------

          unit_func("face()");

          int ix1,iy1,iz1;
          int ix2,iy2,iz2;

          face_1->get_face(&ix1,&iy1,&iz1);
          face_2->get_face(&ix2,&iy2,&iz2);

          unit_assert((ix1 == -ix2));
          unit_assert((iy1 == -iy2));
          unit_assert((iz1 == -iz2));

          unit_assert(std::abs(ix1) + std::abs(ix2) == 2*(iaxis==0)?1:0);
          unit_assert(std::abs(iy1) + std::abs(iy2) == 2*(iaxis==1)?1:0);
          unit_assert(std::abs(iz1) + std::abs(iz2) == 2*(iaxis==2)?1:0);

          //--------------------------------------------------

          delete face_1;
          delete face_2;

        } // ilevel
      } // iface
    } // iaxis
  } // rank

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

