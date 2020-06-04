// See LICENSE_CELLO file for license and copyright information

/// @file     test_FaceFluxes.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    Test program for the FaceFluxes class

#include "main.hpp"
#include "test.hpp"

#include "data.hpp"

#include "test_setup_face.hpp"

#define DEBUG_FACE_FLUXES

double init_1(int ix, int iy, int iz, int mx, int my, int mz)
{
  return 7.0+3.0*(ix+mx*(iy+my*iz));
  //    return 1.0;
}

double init_2(int ix, int iy, int iz, int mx, int my, int mz)
{
  return 1.0+2.0*(ix+mx*(iy+my*iz));
  //  return 2.0;
}

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("FaceFluxes");

  for (int index_test=0; index_test<test::num_face_tests; index_test++) {

    CkPrintf ("Running test %d\n",index_test);

    const int test = test::face_test[index_test].test;
    const int level_1 = test::face_test[index_test].levels_1;
    const int level_2 = test::face_test[index_test].levels_2;
    const int rank = test::face_test[index_test].rank;
    const int ax = test::face_test[index_test].array[0];
    const int ay = test::face_test[index_test].array[1];
    const int az = test::face_test[index_test].array[2];
    const int px = test::face_test[index_test].periodic[0];
    const int py = test::face_test[index_test].periodic[1];
    const int pz = test::face_test[index_test].periodic[2];
    const int NX = test::face_test[index_test].NX;
    const int NY = test::face_test[index_test].NY;
    const int NZ = test::face_test[index_test].NZ;
    const int gx = test::face_test[index_test].ghost[0];
    const int gy = test::face_test[index_test].ghost[1];
    const int gz = test::face_test[index_test].ghost[2];
    const int cx = test::face_test[index_test].cx;
    const int cy = test::face_test[index_test].cy;
    const int cz = test::face_test[index_test].cz;
    const double hx1 = test::face_test[index_test].hx;
    const double hy1 = test::face_test[index_test].hy;
    const double hz1 = test::face_test[index_test].hz;
    const double dt1 =  test::face_test[index_test].dt;
    double hx2 = hx1;
    double hy2 = hy1;
    double hz2 = hz1;
    double dt2 = dt1;
    if (level_1 < level_2) { hx2*=0.5; hy2*=0.5; hz2*=0.5; dt2*=0.5;}
    if (level_1 > level_2) { hx2*=2.0; hy2*=2.0; hz2*=2.0; dt2*=2.0;}
    const int face = test::face_test[index_test].face;
    const int axis = test::face_test[index_test].axis;

    const int fx = 0;
    const int fy = 0;
    const int fz = 0;
    
    Index index_1 = test::create_index
      ( test::face_test[index_test].levels_1,
        test::face_test[index_test].array_1,
        test::face_test[index_test].tree_1);
    Index index_2 = test::create_index
      ( test::face_test[index_test].levels_2,
        test::face_test[index_test].array_2,
        test::face_test[index_test].tree_2);
    
    Face * face_1 = new Face (index_1,index_2,rank,ax,ay,az,px,py,pz,fx,fy,fz);
    Face * face_2 = new Face (index_2,index_1,rank,ax,ay,az,px,py,pz,fx,fy,fz);
    face_1->set_normal(axis,face);
    face_2->set_normal(axis,-face);

  
    int index_field = 3;
    auto face_fluxes_1 = new FaceFluxes
      (*face_1,index_field, NX,NY,NZ, hx1,hy1,hz1,dt1);
    auto face_fluxes_2 = new FaceFluxes
      (*face_2,index_field, NX,NY,NZ, hx2,hy2,hz2,dt2);

    unit_func ("test number");
    unit_assert(index_test == test);

    unit_assert (face_fluxes_1 != NULL);
    unit_assert (face_fluxes_2 != NULL);

    unit_func ("get_dimensions()");
    
    int mx,my,mz;
    face_fluxes_1->get_dimensions (&mx,&my,&mz);

    CkPrintf ("DEBUG_GET_DIMENSIONS (mx,my,mz) %d %d %d\n",mx,my,mz);
    CkPrintf ("DEBUG_GET_DIMENSIONS (NX,NY,NZ) %d %d %d\n",NX,NY,NZ);
    CkPrintf ("DEBUG_GET_DIMENSIONS axis %d\n",axis);

    unit_assert (mx==((axis == 0) ? 1 : NX));
    unit_assert (my==((axis == 1) ? 1 : NY));
    unit_assert (mz==((axis == 2) ? 1 : NZ));

    unit_func ("get_limits()");
    
    int ixl,iyl,izl;
    int ixu,iyu,izu;

    if (level_1 <= level_2) {
      face_fluxes_1 -> get_limits
        (&ixl,&ixu,&iyl,&iyu,&izl,&izu);
      CkPrintf ("DEBUG_LIMITS %d %d %d  %d %d %d\n",
                (ixu-ixl),(iyu-iyl),(izu-izl) ,mx,my,mz);

      // unit_assert ( (ixu-ixl)*(iyu-iyl)*(izu-izl) == mx*my*mz);
      unit_assert ( unit_incomplete );
    }
                                   
    unit_func ("set_ghost()");

    face_fluxes_1->set_ghost(gx,gy,gz);
    face_fluxes_2->set_ghost(gx,gy,gz);

    face_fluxes_1->get_dimensions (&mx,&my,&mz);

    unit_assert (mx==((axis == 0) ? 1 : NX + 2*gx));
    unit_assert (my==((axis == 1) ? 1 : (rank>=2)?NY + 2*gy : 1));
    unit_assert (mz==((axis == 2) ? 1 : (rank>=3)?NZ + 2*gz : 1));
    CkPrintf ("DEBUG_SET_GHOST %d %d %d %d  (%d %d %d) + 2*(%d %d %d)\n",
              __LINE__,mx,my,mz,NX,NY,NZ,gx,gy,gz);

    unit_func("set_centering()");

    face_fluxes_1->set_centering(cx,cy,cz);
    face_fluxes_2->set_centering(cx,cy,cz);

    face_fluxes_1->get_dimensions (&mx,&my,&mz);
    CkPrintf ("DEBUG_COARSEN %d %d %d %d\n",__LINE__,mx,my,mz);
    unit_assert (mx==((axis == 0) ? 1 : NX + 2*gx + cx));
    unit_assert (my==((axis == 1) ? 1 : (rank>=2)?(NY + 2*gy + cy):1));
    unit_assert (mz==((axis == 2) ? 1 : (rank>=3)?(NZ + 2*gz + cz):1));

    

    unit_func ("allocate()");

    face_fluxes_1->allocate();
    face_fluxes_2->allocate();
    

    unit_func("set_fluxes()"); 

    std::vector<double> array_1;
    long int sum_1 = 0;
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
    face_fluxes_1->set_fluxes(array_1,my*mz,mz,1);
    
    std::vector<double> array_2;
    long int sum_2 = 0;
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
    face_fluxes_2->set_fluxes(array_2,my*mz,mz,1);


    int dx1,dy1,dz1;
    auto fluxes_1 = face_fluxes_1->get_fluxes(&dx1,&dy1,&dz1);
    int dx2,dy2,dz2;
    auto fluxes_2 = face_fluxes_2->get_fluxes(&dx2,&dy2,&dz2);
    unit_assert (fluxes_1.size() == mx*my*mz);
    unit_assert (fluxes_2.size() == mx*my*mz);

    unit_func ("get_fluxes()");
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

    // void set_centering(int cx, int cy, int cz)

    unit_func("face()"); 

    unit_assert (face_fluxes_1->face() == *face_1);
    unit_assert (face_fluxes_2->face() == *face_2);

    unit_func("get_element_size()"); 

    double hx1o,hy1o,hz1o;
    face_fluxes_1->get_element_size (&hx1o, &hy1o, &hz1o);

    unit_assert (hx1 == hx1o);
    unit_assert (hy1 == hy1o);
    unit_assert (hz1 == hz1o);

    double hx2o,hy2o,hz2o;
    face_fluxes_2->get_element_size (&hx2o, &hy2o, &hz2o);

    unit_assert (hx2 == hx2o);
    unit_assert (hy2 == hy2o);
    unit_assert (hz2 == hz2o);

    CkPrintf ("DEBUG_COARSEN hx hy hz 1  %g %g %g\n",hx1o,hy1o,hz1o);
    CkPrintf ("DEBUG_COARSEN hx hy hz 2  %g %g %g\n",hx2o,hy2o,hz2o);
    unit_func("time_step()"); 

    unit_assert (dt1 == face_fluxes_1->time_step());
    unit_assert (dt2 == face_fluxes_2->time_step());
    unit_assert (((level_1 == level_2) && (dt1 == dt2)) ||
                 ((level_1 != level_2) && (dt1 != dt2)));
    unit_assert (((level_1 == level_2) && (hx1 == hx2)) ||
                 ((level_1 != level_2) && (hx1 != hx2)));
    unit_assert (rank < 2 ||
                 (((level_1 == level_2) && (hy1 == hy2)) ||
                  ((level_1 != level_2) && (hy1 != hy2))));
    unit_assert (rank < 3 ||
                 (((level_1 == level_2) && (hz1 == hz2)) ||
                  ((level_1 != level_2) && (hz1 != hz2))));

    unit_func("float ratio_cell_width()");
  
    unit_assert (ratio_cell_width(*face_fluxes_1,*face_fluxes_2) == hx1/hx2);
  
    unit_func("float ratio_time_step()");
  
    unit_assert (ratio_time_step(*face_fluxes_1,*face_fluxes_2) == dt1/dt2);

    // --------------------------------------------------
    
    //  void coarsen ()

    unit_func("coarsen()");   

    // weighting for summing: # fine face cells in a coarse face cell per rank
    
    const long w3[] = {0,1,2,4};
    
    if (level_1 > level_2) {

      face_fluxes_1->get_dimensions (&mx,&my,&mz);
      CkPrintf ("DEBUG coarsening face_fluxes_1 before %d %d %d\n",mx,my,mz);

      face_fluxes_1->coarsen();

      face_fluxes_1->get_dimensions (&mx,&my,&mz);
      CkPrintf ("DEBUG coarsening face_fluxes_1 after  %d %d %d\n",mx,my,mz);

      unit_assert (mx==((axis == 0) ? 1 : NX/2 + 2*gx + cx));
      unit_assert (my==((axis == 1) ? 1 : NY/2 + 2*gy + cy));
      unit_assert (mz==((axis == 2) ? 1 : NZ/2 + 2*gz + cz));

      int dxc,dyc,dzc;
      auto fluxes_coarse_1 = face_fluxes_1->get_fluxes(&dxc,&dyc,&dzc);
      long sum_1_coarse = 0;
      for (int iz=0; iz<mz; iz++) {
        for (int iy=0; iy<my; iy++) {
          for (int ix=0; ix<mx; ix++) {
            int i=ix*dxc+iy*dyc+iz*dzc;
            sum_1_coarse += fluxes_coarse_1[i];
          }
        }
      }

      unit_assert (sum_1*w3[rank] == sum_1_coarse);
      
      
    } else if (level_2 > level_1) {

      face_fluxes_2->get_dimensions (&mx,&my,&mz);
      CkPrintf ("DEBUG coarsening face_fluxes_2 before %d %d %d\n",mx,my,mz);

      face_fluxes_2->coarsen();

      face_fluxes_2->get_dimensions (&mx,&my,&mz);
      CkPrintf ("DEBUG coarsening face_fluxes_2 after  %d %d %d\n",mx,my,mz);

      unit_assert (mx==((axis == 0) ? 1 : NX/2 + 2*gx + cx));
      unit_assert (my==((axis == 1) ? 1 : NY/2 + 2*gy + cy));
      unit_assert (mz==((axis == 2) ? 1 : NZ/2 + 2*gz + cz));

      int dxc,dyc,dzc;
      auto fluxes_coarse_2 = face_fluxes_2->get_fluxes(&dxc,&dyc,&dzc);
      long sum_2_coarse = 0;
      for (int iz=0; iz<mz; iz++) {
        for (int iy=0; iy<my; iy++) {
          for (int ix=0; ix<mx; ix++) {
            int i=ix*dxc+iy*dyc+iz*dzc;
            sum_2_coarse += fluxes_coarse_2[i];
          }
        }
      }

      unit_assert (sum_2*w3[rank] == sum_2_coarse);
      
    }

    unit_assert (ratio_cell_width(*face_fluxes_1,*face_fluxes_2) == 1.0);

    //--------------------------------------------------

    //  FaceFluxes & operator *= (double)

    {
      unit_func("operator  *=()");

      auto & fluxes_1 = face_fluxes_1->get_fluxes(&dx1,&dy1,&dz1);
      int mx1,my1,mz1;
      face_fluxes_1->get_dimensions (&mx1,&my1,&mz1);
  
      long sum_1_pre = 0;
      for (int iz=0; iz<mz1; iz++) {
        for (int iy=0; iy<my1; iy++) {
          for (int ix=0; ix<mx1; ix++) {
            int i=ix*dx1+iy*dy1+iz*dz1;
            sum_1_pre += fluxes_1[i];
          }
        }
      }

      (*face_fluxes_1) *= (w3[rank]);

      long sum_1_post = 0;
      for (int iz=0; iz<mz1; iz++) {
        for (int iy=0; iy<my1; iy++) {
          for (int ix=0; ix<mx1; ix++) {
            int i=ix*dx1+iy*dy1+iz*dz1;
            sum_1_post += fluxes_1[i];
          }
        }
      }

      unit_assert (sum_1_post == sum_1_pre * w3[rank]);
    }
    
    //--------------------------------------------------
    
    // FaceFluxes & operator += ()
    
    unit_func("operator += ()");

    {

      auto & fluxes_1 = face_fluxes_1->get_fluxes(&dx1,&dy1,&dz1);
      auto & fluxes_2 = face_fluxes_2->get_fluxes(&dx2,&dy2,&dz2);

      int mx1,my1,mz1;
      int ox1,oy1,oz1;
      face_fluxes_1->get_dimensions (&mx1,&my1,&mz1);
      face_1->get_offset (&ox1,&oy1,&oz1,NX,NY,NZ);
      int mx2,my2,mz2;
      int ox2,oy2,oz2;
      face_fluxes_2->get_dimensions (&mx2,&my2,&mz2);
      face_2->get_offset (&ox2,&oy2,&oz2,NX,NY,NZ);

#ifdef DEBUG_FACE_FLUXES
      CkPrintf ("DEBUG_FACE_FLUXES dimensions %d %d %d  %d %d %d\n",
                mx1,my1,mz1,mx2,my2,mz2);
      CkPrintf ("DEBUG_FACE_FLUXES offset %d %d %d  %d %d %d\n",
                ox1,oy1,oz1,ox2,oy2,oz2);
#endif      
      long sum_1_pre = 0;
      for (int iz=0; iz<mz1; iz++) {
        for (int iy=0; iy<my1; iy++) {
          for (int ix=0; ix<mx1; ix++) {
            int i=ix*dx1+iy*dy1+iz*dz1;
            sum_1_pre += fluxes_1[i];
          }
        }
      }
      long sum_2_pre = 0;
      for (int iz=0; iz<mz2; iz++) {
        for (int iy=0; iy<my2; iy++) {
          for (int ix=0; ix<mx2; ix++) {
            int i=ix*dx2+iy*dy2+iz*dz2;
            sum_2_pre += fluxes_2[i];
          }
        }
      }

      (*face_fluxes_1) += (*face_fluxes_2);

      long sum_1_post = 0;
      for (int iz=0; iz<mz1; iz++) {
        for (int iy=0; iy<my1; iy++) {
          for (int ix=0; ix<mx1; ix++) {
            int i=ix*dx1+iy*dy1+iz*dz1;
            sum_1_post += fluxes_1[i];
          }
        }
      }

      // unit_assert (sum_1_post == (sum_1_pre + sum_2_pre));
      unit_assert ( unit_incomplete );
#ifdef DEBUG_FACE_FLUXES
      CkPrintf ("sum_1 %ld =?= %ld (sum_1_pre %ld + sum_2_pre %ld)\n",
                sum_1_post,(sum_1_pre+sum_2_pre),sum_1_pre,sum_2_pre);
#endif        
      //      unit_assert (unit_incomplete);
    }

    //--------------------------------------------------

    //  friend FaceFluxes operator - (FaceFluxes ff_1, FaceFluxes ff_2)

    unit_func("FaceFluxes operator -()");
  
    unit_assert (unit_incomplete);

    //--------------------------------------------------

    // deallocate()
    
    unit_func("deallocate()");

    face_fluxes_1->deallocate();
    face_fluxes_2->deallocate();

    unit_assert (face_fluxes_1->get_fluxes(0,0,0).size() == 0);
    unit_assert (face_fluxes_2->get_fluxes(0,0,0).size() == 0);

    delete face_fluxes_1;
    delete face_fluxes_2;
  }
  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

