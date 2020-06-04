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

  const int ax=4, ay=4, az=4; // array size
  const int tx=3, ty=1, tz=2; // tree in array

  int index_field = 3;
  const int nx = 16;
  const int ny = 16;
  const int nz = 16;
  const double dt = 0.125;
  int cx=1, cy=0, cz=1;
  FaceFluxes * ff27_blk[3][3][3];
  FaceFluxes * ff27_nbr[3][3][3];


  unit_func ("deallocate_fluxes()");
  FluxData flux_data;
  for (int iz=-1; iz<=1; iz++) {
    for (int iy=-1; iy<=1; iy++) {
      for (int ix=-1; ix<=1; ix++) {
        ff27_blk[ix+1][iy+1][iz+1] = nullptr;
        ff27_nbr[ix+1][iy+1][iz+1] = nullptr;
      }
    }
  }
  // create face fluxes
  const int level = 3;
  const Face f011(0,1,1,1,0,0,cx,cy,cz);
  const Face f211(2,1,1,1,0,0,cx,cy,cz);
  const Face f101(1,0,1,0,1,0,cx,cy,cz);
  const Face f121(1,2,1,0,1,0,cx,cy,cz);
  const Face f110(1,1,0,0,0,1,cx,cy,cz);
  const Face f112(1,1,2,0,0,1,cx,cy,cz);

  unit_func("fluxes_blk");

  ff27_blk[0][1][1] = new FaceFluxes (f011,index_field,nx,ny,nz,level,dt);
  ff27_blk[2][1][1] = new FaceFluxes (f211,index_field,nx,ny,nz,level,dt);
  ff27_blk[1][0][1] = new FaceFluxes (f101,index_field,nx,ny,nz,level,dt);
  ff27_blk[1][2][1] = new FaceFluxes (f121,index_field,nx,ny,nz,level,dt);
  ff27_blk[1][1][0] = new FaceFluxes (f110,index_field,nx,ny,nz,level,dt);
  ff27_blk[1][1][2] = new FaceFluxes (f112,index_field,nx,ny,nz,level,dt);

  ff27_nbr[0][1][1] = new FaceFluxes (f011,index_field,nx,ny,nz,level,dt);
  ff27_nbr[2][1][1] = new FaceFluxes (f211,index_field,nx,ny,nz,level,dt);
  ff27_nbr[1][0][1] = new FaceFluxes (f101,index_field,nx,ny,nz,level,dt);
  ff27_nbr[1][2][1] = new FaceFluxes (f121,index_field,nx,ny,nz,level,dt);
  ff27_nbr[1][1][0] = new FaceFluxes (f110,index_field,nx,ny,nz,level,dt);
  ff27_nbr[1][1][2] = new FaceFluxes (f112,index_field,nx,ny,nz,level,dt);

  unit_func ("FluxData::fluxes_block()");

  unit_assert(flux_data.fluxes_block(-1,0,0) == nullptr);
  unit_assert(flux_data.fluxes_block(+1,0,0) == nullptr);
  unit_assert(flux_data.fluxes_block(0,-1,0) == nullptr);
  unit_assert(flux_data.fluxes_block(0,+1,0) == nullptr);
  unit_assert(flux_data.fluxes_block(0,0,-1) == nullptr);
  unit_assert(flux_data.fluxes_block(0,0,+1) == nullptr);
  
  unit_assert(flux_data.fluxes_neighbor(-1,0,0) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(+1,0,0) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(0,-1,0) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(0,+1,0) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(0,0,-1) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(0,0,+1) == nullptr);

  unit_func ("FluxData::set_fluxes_block()");
  
  flux_data.set_fluxes_block(ff27_blk[0][1][1],-1,0,0);
  flux_data.set_fluxes_block(ff27_blk[2][1][1],+1,0,0);
  flux_data.set_fluxes_block(ff27_blk[1][0][1],0,-1,0);
  flux_data.set_fluxes_block(ff27_blk[1][2][1],0,+1,0);
  flux_data.set_fluxes_block(ff27_blk[1][1][0],0,0,-1);
  flux_data.set_fluxes_block(ff27_blk[1][1][2],0,0,+1);

  unit_assert(flux_data.fluxes_block(-1,0,0) == ff27_blk[0][1][1]);
  unit_assert(flux_data.fluxes_block(+1,0,0) == ff27_blk[2][1][1]);
  unit_assert(flux_data.fluxes_block(0,-1,0) == ff27_blk[1][0][1]);
  unit_assert(flux_data.fluxes_block(0,+1,0) == ff27_blk[1][2][1]);
  unit_assert(flux_data.fluxes_block(0,0,-1) == ff27_blk[1][1][0]);
  unit_assert(flux_data.fluxes_block(0,0,+1) == ff27_blk[1][1][2]);
  
  unit_func ("FluxData::set_fluxes_neighbor()");
  
  flux_data.set_fluxes_neighbor(ff27_nbr[0][1][1],-1,0,0);
  flux_data.set_fluxes_neighbor(ff27_nbr[2][1][1],+1,0,0);
  flux_data.set_fluxes_neighbor(ff27_nbr[1][0][1],0,-1,0);
  flux_data.set_fluxes_neighbor(ff27_nbr[1][2][1],0,+1,0);
  flux_data.set_fluxes_neighbor(ff27_nbr[1][1][0],0,0,-1);
  flux_data.set_fluxes_neighbor(ff27_nbr[1][1][2],0,0,+1);

  unit_assert(flux_data.fluxes_neighbor(-1,0,0) == ff27_nbr[0][1][1]);
  unit_assert(flux_data.fluxes_neighbor(+1,0,0) == ff27_nbr[2][1][1]);
  unit_assert(flux_data.fluxes_neighbor(0,-1,0) == ff27_nbr[1][0][1]);
  unit_assert(flux_data.fluxes_neighbor(0,+1,0) == ff27_nbr[1][2][1]);
  unit_assert(flux_data.fluxes_neighbor(0,0,-1) == ff27_nbr[1][1][0]);
  unit_assert(flux_data.fluxes_neighbor(0,0,+1) == ff27_nbr[1][1][2]);

  unit_func ("FluxData::delete_fluxes_block()");

  flux_data.delete_fluxes_block (-1,0,0);
  flux_data.delete_fluxes_block (+1,0,0);
  flux_data.delete_fluxes_block (0,-1,0);
  flux_data.delete_fluxes_block (0,+1,0);
  flux_data.delete_fluxes_block (0,0,-1);
  flux_data.delete_fluxes_block (0,0,+1);
  
  unit_assert(flux_data.fluxes_block(-1,0,0) == nullptr);
  unit_assert(flux_data.fluxes_block(+1,0,0) == nullptr);
  unit_assert(flux_data.fluxes_block(0,-1,0) == nullptr);
  unit_assert(flux_data.fluxes_block(0,+1,0) == nullptr);
  unit_assert(flux_data.fluxes_block(0,0,-1) == nullptr);
  unit_assert(flux_data.fluxes_block(0,0,+1) == nullptr);
  
  unit_func ("FluxData::delete_fluxes_neighbor()");

  flux_data.delete_fluxes_neighbor (-1,0,0);
  flux_data.delete_fluxes_neighbor (+1,0,0);
  flux_data.delete_fluxes_neighbor (0,-1,0);
  flux_data.delete_fluxes_neighbor (0,+1,0);
  flux_data.delete_fluxes_neighbor (0,0,-1);
  flux_data.delete_fluxes_neighbor (0,0,+1);
  
  unit_assert(flux_data.fluxes_neighbor(-1,0,0) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(+1,0,0) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(0,-1,0) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(0,+1,0) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(0,0,-1) == nullptr);
  unit_assert(flux_data.fluxes_neighbor(0,0,+1) == nullptr);

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

