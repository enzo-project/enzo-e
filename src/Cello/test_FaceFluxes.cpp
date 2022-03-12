// See LICENSE_CELLO file for license and copyright information

/// @file     test_FaceFluxes.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    Test program for the FaceFluxes class

#include "main.hpp"
#include "test.hpp"

#include "data.hpp"

cello_float init_1(int ix, int iy, int iz, int mx, int my, int mz)
{
  return 7.0+3.0*(ix+mx*(iy+my*iz));
}

cello_float init_2(int ix, int iy, int iz, int mx, int my, int mz)
{
  return 1.0+2.0*(ix+mx*(iy+my*iz));
}

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("FaceFluxes");

  const int num_level = 3;
  const int level_1[] = {1,2,3};
  const int level_2[] = {1,3,2};

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

  const int rvol[4] = {1, 2, 4, 8};
  
  for (int rank = 2; rank <=3; ++rank) {
    for (int iaxis=0; iaxis<rank; iaxis++) {
      for (int iface=0; iface<2; ++iface) {
        for (int ilevel=0; ilevel<num_level; ++ilevel) {
        
          const int L_1 = level_1[ilevel];
          const int L_2 = level_2[ilevel];
          int n3[3];
          n3[0] = test.size[0];
          n3[1] = (rank >= 2) ? test.size[1]:1;
          n3[2] = (rank >= 3) ? test.size[2]:1;
          int c3[3];
          c3[0] = test.centered[0];
          c3[1] = (rank >= 2) ? test.centered[1] : 0;
          c3[2] = (rank >= 3) ? test.centered[2] : 0;
          double h3_1[3];
          h3_1[0] = test.cell_width[0];
          h3_1[1] = (rank >= 2) ? test.cell_width[1] : 0;
          h3_1[2] = (rank >= 3) ? test.cell_width[2] : 0;
          double h3_2[3] = {h3_1[0],h3_1[1],h3_1[2]};
          const int cx = test.child[0];
          const int cy = (rank >= 2) ? test.child[1] : 0;
          const int cz = (rank >= 3) ? test.child[2] : 0;

          if (L_1 < L_2) { h3_2[0]*=0.5; h3_2[1]*=0.5; h3_2[2]*=0.5;}
          if (L_1 > L_2) { h3_2[0]*=2.0; h3_2[1]*=2.0; h3_2[2]*=2.0;}
          int fx = iaxis==0 ? iface*2-1 : 0;
          int fy = iaxis==1 ? iface*2-1 : 0;
          int fz = iaxis==2 ? iface*2-1 : 0;

          Face face_1( fx, fy, fz,iaxis,iface);
          Face face_2(-fx,-fy,-fz,iaxis,iface);
  
          int index_field = 3;
          auto face_fluxes_1 = new FaceFluxes
            (face_1,index_field, n3[0],n3[1],n3[2], c3[0],c3[1],c3[2]);
          auto face_fluxes_2 = new FaceFluxes
            (face_2,index_field, n3[0],n3[1],n3[2], c3[0],c3[1],c3[2]);

          unit_assert (face_fluxes_1 != NULL);
          unit_assert (face_fluxes_2 != NULL);

          unit_func ("get_size()");
    
          int mx,my,mz;
          face_fluxes_1->get_size (&mx,&my,&mz);

          unit_assert (mx == ((fx!=0) ? 1 : n3[0]));
          if (rank >= 2) unit_assert (my == ((fy!=0) ? 1 : n3[1]));
          if (rank >= 3) unit_assert (mz == ((fz!=0) ? 1 : n3[2]));

          int nx1,ny1,nz1;

          face_fluxes_1 -> get_size  (&nx1, &ny1, &nz1);
    
    
          if (rank == 3) {
            // 3D
            int a1 = (iaxis+1) % 3;
            int a2 = (iaxis+2) % 3;
            unit_assert ( nx1*ny1*nz1 == n3[a1]*n3[a2]);
          } else if (rank == 2) {
            // 2D
            int a1 = (iaxis+1) % 2;
            unit_assert ( nx1*ny1*nz1 == n3[a1] );
          }

          //--------------------------------------------------

          int nx2,ny2,nz2;
          face_fluxes_2 -> get_size  (&nx2, &ny2, &nz2);

          if (rank >= 3) {
            int a1 = (iaxis+1) % 3;
            int a2 = (iaxis+2) % 3;
            unit_assert ( nx2*ny2*nz2 == n3[a1]*n3[a2]);
          } else if (rank >= 2) {
            int a1 = (iaxis+1) % 2;
            unit_assert ( nx2*ny2*nz2 == n3[a1]);
          }

          face_fluxes_1->get_size(&mx,&my,&mz);

          unit_assert (mx == ((fx!=0) ? 1 : n3[0]));
          unit_assert (my == ((fy!=0) ? 1 : (rank>=2)?n3[1] : 1));
          unit_assert (mz == ((fz!=0) ? 1 : (rank>=3)?n3[2] : 1));

          face_fluxes_1->get_size(&mx,&my,&mz);

          unit_assert (mx == ((fx!=0) ? 1 : n3[0] + c3[0]));
          unit_assert (my == ((fy!=0) ? 1 : (rank>=2)?(n3[1] + c3[1]):1));
          unit_assert (mz == ((fz!=0) ? 1 : (rank>=3)?(n3[2] + c3[2]):1));


          unit_func ("allocate()");

          unit_assert(face_fluxes_1->flux_array() == nullptr);
          unit_assert(face_fluxes_2->flux_array() == nullptr);
          face_fluxes_1->allocate_storage();
          face_fluxes_2->allocate_storage();
          unit_assert(face_fluxes_1->get_size() == mx*my*mz);
          unit_assert(face_fluxes_2->get_size() == mx*my*mz);

          unit_func("set_flux_array()"); 

          std::vector<cello_float> array_1;
          double sum_1 = 0;
          array_1.resize(mx*my*mz);
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                int i=iz+mz*(iy+my*ix);
                array_1[i] = init_1(ix,iy,iz,mx,my,mz);
                sum_1 += array_1[i];
              }
            }
          }
          face_fluxes_1->set_flux_array(array_1,my*mz,mz,1);
    
          std::vector<cello_float> array_2;
          double sum_2 = 0;
          array_2.resize(mx*my*mz);
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                int i=iz+mz*(iy+my*ix);
                array_2[i] = init_2(ix,iy,iz,mx,my,mz);
                sum_2 += array_2[i];
              }
            }
          }
          face_fluxes_2->set_flux_array(array_2,my*mz,mz,1);


          int dx1,dy1,dz1;
          cello_float * fluxes_1 = face_fluxes_1->flux_array(&dx1,&dy1,&dz1);
          int dx2,dy2,dz2;
          cello_float * fluxes_2 = face_fluxes_2->flux_array(&dx2,&dy2,&dz2);

          unit_func ("flux_array()");
          bool match_1 = true, match_2=true;
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                int i_f1=ix*dx1 + iy*dy1 + iz*dz1;
                int i_f2=ix*dx2 + iy*dy2 + iz*dz2;
                int i_a1=iz+mz*(iy+my*ix);
                int i_a2=iz+mz*(iy+my*ix);
                if ( fluxes_1[i_f1] != array_1[i_a1]) {
                  match_1 = false;
                }
                if ( fluxes_2[i_f2] != array_2[i_a2]) {
                  match_2 = false;
                }
              }
            }
          }
          unit_assert (match_1);
          unit_assert (match_2);

          // --------------------------------------------------
    
          //  void coarsen ()

          unit_func("coarsen()");   

          if (L_1 > L_2) {

            face_fluxes_1->coarsen(cx,cy,cz,rank);

            face_fluxes_1->get_size (&mx,&my,&mz);

            int dxc,dyc,dzc;
            auto fluxes_coarse_1 = face_fluxes_1->flux_array(&dxc,&dyc,&dzc);
            double sum_1_coarse = 0;
            for (int iz=0; iz<mz; iz++) {
              for (int iy=0; iy<my; iy++) {
                for (int ix=0; ix<mx; ix++) {
                  int i=ix*dxc+iy*dyc+iz*dzc;
                  sum_1_coarse += fluxes_coarse_1[i];
                }
              }
            }

            unit_assert (sum_1 == sum_1_coarse*rvol[rank]);

          } else if (L_2 > L_1) {

            face_fluxes_2->get_size(&mx,&my,&mz);
            face_fluxes_2->coarsen(cx,cy,cz,rank);

            int dxc,dyc,dzc;
            auto fluxes_coarse_2 = face_fluxes_2->flux_array(&dxc,&dyc,&dzc);
            double sum_2_coarse = 0;
            for (int iz=0; iz<mz; iz++) {
              for (int iy=0; iy<my; iy++) {
                for (int ix=0; ix<mx; ix++) {
                  int i=ix*dxc+iy*dyc+iz*dzc;
                  sum_2_coarse += fluxes_coarse_2[i];
                }
              }
            }

            unit_assert (sum_2 == sum_2_coarse*rvol[rank]);
      
          }

          //  FaceFluxes & operator *= (double)

          {
            unit_func("operator  *=()");

            cello_float * fluxes_1 = face_fluxes_1->flux_array(&dx1,&dy1,&dz1);
            int mx1,my1,mz1;
            face_fluxes_1->get_size (&mx1,&my1,&mz1);
  
            double sum_1_pre = 0;
            for (int iz=0; iz<mz1; iz++) {
              for (int iy=0; iy<my1; iy++) {
                for (int ix=0; ix<mx1; ix++) {
                  int i=ix*dx1+iy*dy1+iz*dz1;
                  sum_1_pre += fluxes_1[i];
                }
              }
            }

            (*face_fluxes_1) *= (7.25);

            double sum_1_post = 0;
            for (int iz=0; iz<mz1; iz++) {
              for (int iy=0; iy<my1; iy++) {
                for (int ix=0; ix<mx1; ix++) {
                  int i=ix*dx1+iy*dy1+iz*dz1;
                  sum_1_post += fluxes_1[i];
                }
              }
            }

            unit_assert (sum_1_post == sum_1_pre * 7.25);
          }
    
          //--------------------------------------------------
    
          // FaceFluxes & operator += ()
    
          unit_func("operator += ()");

          {

            cello_float * fluxes_1 = face_fluxes_1->flux_array(&dx1,&dy1,&dz1);
            cello_float * fluxes_2 = face_fluxes_2->flux_array(&dx2,&dy2,&dz2);

            int mx1,my1,mz1;
            int mx2,my2,mz2;
            face_fluxes_1->get_size (&mx1,&my1,&mz1);
            face_fluxes_2->get_size (&mx2,&my2,&mz2);

            unit_func("operator +=");      

            double sum_1_pre = 0;
            for (int iz=0; iz<mz1; iz++) {
              for (int iy=0; iy<my1; iy++) {
                for (int ix=0; ix<mx1; ix++) {
                  int i=ix*dx1+iy*dy1+iz*dz1;
                  sum_1_pre += fluxes_1[i];
                }
              }
            }
            double sum_2_pre = 0;
            for (int iz=0; iz<mz2; iz++) {
              for (int iy=0; iy<my2; iy++) {
                for (int ix=0; ix<mx2; ix++) {
                  int i=ix*dx2+iy*dy2+iz*dz2;
                  sum_2_pre += fluxes_2[i];
                }
              }
            }

            double sum_1_post = 0;
            double sum_2_post = 0;

            if (L_1 <= L_2) {

              face_fluxes_1 -> accumulate(*face_fluxes_2,cx,cy,cz, rank);

              for (int iz=0; iz<mz1; iz++) {
                for (int iy=0; iy<my1; iy++) {
                  for (int ix=0; ix<mx1; ix++) {
                    int i=ix*dx1+iy*dy1+iz*dz1;
                    sum_1_post += fluxes_1[i];
                  }
                }
              }
            }
            if (L_1 >= L_2) {

              face_fluxes_2->accumulate(*face_fluxes_1,cx,cy,cz,rank);

              for (int iz=0; iz<mz2; iz++) {
                for (int iy=0; iy<my2; iy++) {
                  for (int ix=0; ix<mx2; ix++) {
                    int i=ix*dx2+iy*dy2+iz*dz2;
                    sum_2_post += fluxes_2[i];
                  }
                }
              }
            }

            if (L_1 < L_2) {
              unit_assert (sum_1_post == (sum_1_pre + sum_2_pre));
              unit_assert (sum_2_post == 0);
            } else if (L_1 > L_2) {
              unit_assert (sum_1_post == 0);
              unit_assert (sum_2_post == (sum_1_pre + sum_2_pre));
            } else if (L_1 == L_2) {
              unit_assert (sum_1_post == (sum_1_pre + sum_2_pre));
              unit_assert (sum_2_post == (sum_1_pre + 2*sum_2_pre));
            }
          }

          //--------------------------------------------------
          {
            unit_func("clear()");

            int mx1,my1,mz1;
            int mx2,my2,mz2;
            face_fluxes_1->get_size (&mx1,&my1,&mz1);
            face_fluxes_2->get_size (&mx2,&my2,&mz2);

            {
              double sum_abs = 0.0;
              cello_float * fluxes_1 = face_fluxes_1->flux_array(&dx1,&dy1,&dz1);
              for (int iz=0; iz<mz1; iz++) {
                for (int iy=0; iy<my1; iy++) {
                  for (int ix=0; ix<mx1; ix++) {
                    int i=ix*dx1+iy*dy1+iz*dz1;
                    sum_abs += std::abs(fluxes_1[i]);
                  }
                }
              }
              unit_assert(sum_abs != 0.0);
            }

            face_fluxes_1->clear();
      
            {
              double sum_abs = 0.0;
              cello_float * fluxes_1 = face_fluxes_1->flux_array(&dx1,&dy1,&dz1);
              for (int iz=0; iz<mz1; iz++) {
                for (int iy=0; iy<my1; iy++) {
                  for (int ix=0; ix<mx1; ix++) {
                    int i=ix*dx1+iy*dy1+iz*dz1;
                    sum_abs += std::abs(fluxes_1[i]);
                  }
                }
              }
              unit_assert(sum_abs == 0.0);
            }

            {
              double sum_abs = 0;
              cello_float * fluxes_2 = face_fluxes_2->flux_array(&dx2,&dy2,&dz2);
              for (int iz=0; iz<mz2; iz++) {
                for (int iy=0; iy<my2; iy++) {
                  for (int ix=0; ix<mx2; ix++) {
                    int i=ix*dx2+iy*dy2+iz*dz2;
                    sum_abs += std::abs(fluxes_2[i]);
                  }
                }
              }
              unit_assert(sum_abs != 0.0);
            }
            face_fluxes_2->clear();
            {
              double sum_abs = 0;
              cello_float * fluxes_2 = face_fluxes_2->flux_array(&dx2,&dy2,&dz2);
              for (int iz=0; iz<mz2; iz++) {
                for (int iy=0; iy<my2; iy++) {
                  for (int ix=0; ix<mx2; ix++) {
                    int i=ix*dx2+iy*dy2+iz*dz2;
                    sum_abs += std::abs(fluxes_2[i]);
                  }
                }
              }
              unit_assert(sum_abs == 0.0);
            }
          }
          // deallocate()
    
          unit_func("deallocate()");

          face_fluxes_1->deallocate_storage();
          face_fluxes_2->deallocate_storage();

          unit_assert (face_fluxes_1->flux_array() == nullptr);
          unit_assert (face_fluxes_2->flux_array() == nullptr);

          delete face_fluxes_1;
          delete face_fluxes_2;

        } // ilevel
      } // iface
    } // iaxis
  } // rank

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

