// See LICENSE_CELLO file for license and copyright information

/// @file     test_FaceFluxes.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    Test program for the FaceFluxes class

#include "main.hpp"
#include "test.hpp"

#include "data.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("FaceFluxes");

  Index index_block(3,1,2);
  Index index_neighbor(3,2,2);
  const int ax=4, ay=4, az=4;

  Face face_1(index_block,index_neighbor,3,ax,ay,az);
  Face face_2(index_neighbor,index_block,3,ax,ay,az);
  int index_field = 3;
  const int nx = 16;
  const int ny = 16;
  const int nz = 16;
  const double hx = 0.25;
  const double hy = 1.00;
  const double hz = 0.25;
  const double dt = 0.125;
  auto face_fluxes_1 = new FaceFluxes (face_1,index_field,
                                       nx,ny,nz,
                                       hx,hy,hz,dt);
  auto face_fluxes_2 = new FaceFluxes (face_2,index_field,
                                       nx,ny,nz,
                                       hx,hy,hz,dt);

  unit_assert (face_fluxes_1 != NULL);
  unit_assert (face_fluxes_2 != NULL);

  int mx,my,mz;
  face_fluxes_1->get_dimensions (&mx,&my,&mz);

  unit_assert (mx==nx);
  unit_assert (my==1);
  unit_assert (mz==nz);

  unit_func("set_ghost");

  face_fluxes_1->set_ghost(3,3,3);

  face_fluxes_1->get_dimensions (&mx,&my,&mz);

  unit_assert (mx==nx+2*3);
  unit_assert (my==1);
  unit_assert (mz==nz+2*3);

  unit_func("set_centering");

  face_fluxes_1->set_centering(1,1,0);

  face_fluxes_1->get_dimensions (&mx,&my,&mz);

  unit_assert (mx==nx+2*3+1);
  unit_assert (my==1);
  unit_assert (mz==nz+2*3);

  unit_func ("allocate");
  face_fluxes_1->allocate();

  std::vector<double> array;
  array.resize(mx*my*mz);
  for (int ix=0; ix<mx; ix++) {
    for (int iy=0; iy<my; iy++) {
      for (int iz=0; iz<mz; iz++) {
        int i=iz+mz*(iy+my*ix);
        array[i] = ix+mx*(iy+my*iz);
      }
    }
  }
  face_fluxes_1->set_fluxes(array,my*mz,mz,1);

  for (int ix=0; ix<mx; ix++) {
    for (int iy=0; iy<my; iy++) {
      for (int iz=0; iz<mz; iz++) {
        int i=iz+mz*(iy+my*ix);
        array[i] = ix+mx*(iy+my*iz);
      }
    }
  }

  int dx,dy,dz;
  auto fluxes = face_fluxes_1->get_fluxes(&dx,&dy,&dz);
  unit_assert (fluxes.size() == mx*my*mz);

  unit_func ("get_fluxes");
  bool match = true;
  int jx,jy,jz;
  for (int ix=0; ix<mx; ix++) {
    for (int iy=0; iy<my; iy++) {
      for (int iz=0; iz<mz; iz++) {
        int i_a=iz+mz*(iy+my*ix);
        int i_f=ix*dx + iy*dy + iz*dz;
        if (array[i_a] != fluxes[i_f]) {
          match = false;
          jx=ix;
          jy=iy;
          jz=iz;
        }
      }
    }
  }
  unit_assert (match);

  // void set_centering(int cx, int cy, int cz)

  unit_func("allocate"); 

  //  void allocate ()

  unit_func("deallocate");

  //  void deallocate()

  unit_func("face"); 

  //  Face face () const

  unit_func("get_element_size"); 

  double hxo,hyo,hzo;
  face_fluxes_1->get_element_size (&hxo, &hyo, &hzo);
  unit_assert (hx == hxo);
  unit_assert (hy == hyo);
  unit_assert (hz == hzo);

  unit_func("time_step"); 

  //  double time_step ()

  unit_func("get_dimensions");

  //  void get_dimensions (int *mx, int *my, int *mz) const

  unit_func("set_fluxes"); 

  //  void set_fluxes ( T * array, int dx, int dy, int dz)

  unit_func("float ratio_cell_width");
  
  //  friend float ratio_cell_width
  
  unit_func("float ratio_time_step");
  
  // ratio_time_step (FaceFluxes ff_1, FaceFluxes ff_2)

  unit_func("coarsen");   
  //  void coarsen ()

  unit_func("operator += ()");

  // FaceFluxes & operator += ()

  unit_func("operator  *=");
  
  //  FaceFluxes & operator *= (double)

  unit_func("FaceFluxes operator -");
  
  //  friend FaceFluxes operator - (FaceFluxes ff_1, FaceFluxes ff_2)


  //--------------------------------------------------

  unit_func ("function()");

  //--------------------------------------------------

  delete face_fluxes_1;
  delete face_fluxes_2;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

