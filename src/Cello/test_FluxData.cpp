// See LICENSE_CELLO file for license and copyright information

/// @file     test_FluxData.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-10-15
/// @brief    Test program for the FluxData class

#include "main.hpp"
#include "test.hpp"

#include "data.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("FluxData");

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

  FluxData * flux_data = new FluxData;

  unit_assert (flux_data != nullptr);

  flux_data->insert_fluxes(face_1,face_fluxes_1);
  flux_data->insert_fluxes(face_2,face_fluxes_2);

  //--------------------------------------------------

  unit_func ("function()");

  //--------------------------------------------------

  delete face_fluxes_1;
  delete face_fluxes_2;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

