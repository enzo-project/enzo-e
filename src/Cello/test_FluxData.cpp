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

  // e.g. de vx vy vz
  std::vector<int> field_list = {9,3,4,5};
  std::vector<int> cx_list = {0,1,0,0};
  std::vector<int> cy_list = {0,0,1,0};
  std::vector<int> cz_list = {0,0,0,1};
  const int nf = field_list.size();
  const int n3[3] = {16, 24, 20};
  const double dt = 0.125;
  const int level = 7;

  FluxData flux_data;

  std::vector<double> array_blk[3][2][nf];
  std::vector<double> array_nbr[3][2][nf];

  for (int i_f=0; i_f<nf; i_f++) {

    const int index_field = field_list[i_f];

    for (int axis=0; axis<3; axis++) {

      const int a1 = (axis+1) % 3;
      const int a2 = (axis+2) % 3;

      const int n1=n3[a1];
      const int n2=n3[a2];

      for (int face=0; face<2; face++) {

        FaceFluxes * ff_blk = flux_data.block_fluxes(axis,face,index_field);
        FaceFluxes * ff_nbr = flux_data.neighbor_fluxes(axis,face,index_field);
        
        unit_assert (ff_blk == nullptr);
        unit_assert (ff_nbr == nullptr);

        double value = 7.0 + (axis + 3*(face + 2*(i_f)));

        array_blk[axis][face][i_f].resize(n1*n2);
        array_nbr[axis][face][i_f].resize(n1*n2);

        auto & b = array_blk[axis][face][i_f];
        auto & n = array_nbr[axis][face][i_f];

        for (int i1=0; i1<n1; i1++) {
          for (int i2=0; i2<n2; i2++) {
            int i=i2+n2*(i1);
            b[i] = value + i;
            n[i] = value + i;
          }
        }
      }
    }
  }

  const int n_f = field_list.size();
  
  unit_func ("allocate()");

  flux_data.allocate(n3[0],n3[1],n3[2],level,dt,field_list);
  
  unit_func ("FluxData::set_block_fluxes()");

  for (int i_f=0; i_f<n_f; i_f++) {

    const int index_field = field_list[i_f];

    for (int axis=0; axis<3; axis++) {
      const int n1=n3[(axis+1)%3];
      const int n2=n3[(axis+2)%3];
      for (int face=0; face<2; face++) {

        FaceFluxes * ff_blk = flux_data.block_fluxes(axis,face,index_field);
        FaceFluxes * ff_nbr = flux_data.neighbor_fluxes(axis,face,index_field);

        auto & b = array_blk[axis][face][i_f];
        auto & n = array_nbr[axis][face][i_f];

        unit_assert (ff_blk != nullptr);
        unit_assert (ff_nbr != nullptr);
        
        int mx,my,mz;
        ff_blk->get_dimensions(&mx,&my,&mz);
        unit_func ("get_dimensions()");
        unit_assert(mx*my*mz == n1*n2);

        if (axis==0) {
          ff_blk->set_flux_array (b,0,n3[2],1);
          ff_nbr->set_flux_array (array_nbr[axis][face][i_f],0,n3[2],1);
        } else if (axis==1) {
          ff_blk->set_flux_array (b,1,0,n3[0]);
          ff_nbr->set_flux_array (array_nbr[axis][face][i_f],1,0,n3[0]);
        } else if (axis==2) {
          ff_blk->set_flux_array (b,n3[1],1,0);
          ff_nbr->set_flux_array (array_nbr[axis][face][i_f],n3[1],1,0);
        }

        int dbx,dby,dbz;
        int dnx,dny,dnz;
        std::vector<double> & fluxes_blk = ff_blk->flux_array(&dbx,&dby,&dbz);
        std::vector<double> & fluxes_nbr = ff_nbr->flux_array(&dnx,&dny,&dnz);
        
        int count = 0;
        unit_func ("array_blk()");
        double sum1=0.0,sum2=0.0;
        if (axis==0) {
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                const int i1 = iy;
                const int i2 = iz;
                const int id = ix*dbx+iy*dby+iz*dbz;
                const int ib = i2 + n2*i1;
                sum1+=b[ib];
                sum2+=fluxes_blk[id];
                if (fluxes_blk[id] != b[ib]) count++;
              }
            }
          }
          unit_func("array_blk x-axis");
          unit_assert (count==0);
        } else if (axis==1) {
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                const int i1 = iz;
                const int i2 = ix;
                const int id = ix*dbx+iy*dby+iz*dbz;
                const int ib = i2 + n2*i1;
                sum1+=b[ib];
                sum2+=fluxes_blk[id];
                if (fluxes_blk[id] != b[ib]) count++;
              }
            }
          }
          unit_func("array_blk y-axis");
          unit_assert (count==0);
        } else if (axis==2) {
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                const int i1 = ix;
                const int i2 = iy;
                const int id = ix*dbx+iy*dby+iz*dbz;
                const int ib = i2 + n2*i1;
                sum1+=b[ib];
                sum2+=fluxes_blk[id];
                if (fluxes_blk[id] != b[ib]) count++;
              }
            }
          }
          unit_func("array_blk z-axis");
          unit_assert (count==0);
        }
        unit_func ("array_nbr()");
        count = 0;
        sum1=0.0,sum2=0.0;
        if (axis==0) {
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                const int i1 = iy;
                const int i2 = iz;
                const int id = ix*dbx+iy*dby+iz*dbz;
                const int ib = i2 + n2*i1;
                sum1+=n[ib];
                sum2+=fluxes_nbr[id];
                if (fluxes_nbr[id] != n[ib]) count++;
              }
            }
          }
          unit_func("array_nbr x-axis");
          unit_assert (count==0);
        } else if (axis==1) {
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                const int i1 = iz;
                const int i2 = ix;
                const int id = ix*dbx+iy*dby+iz*dbz;
                const int ib = i2 + n2*i1;
                sum1+=n[ib];
                sum2+=fluxes_nbr[id];
                if (fluxes_nbr[id] != n[ib]) count++;
              }
            }
          }
          unit_func("array_nbr y-axis");
          unit_assert (count==0);
        } else if (axis==2) {
          for (int iz=0; iz<mz; iz++) {
            for (int iy=0; iy<my; iy++) {
              for (int ix=0; ix<mx; ix++) {
                const int i1 = ix;
                const int i2 = iy;
                const int id = ix*dbx+iy*dby+iz*dbz;
                const int ib = i2 + n2*i1;
                sum1+=n[ib];
                sum2+=fluxes_nbr[id];
                if (fluxes_nbr[id] != n[ib]) count++;
              }
            }
          }
          unit_func("array_nbr z-axis");
          unit_assert (count==0);
        }
        
        unit_assert(dbx==dnx && dby==dny && dbz==dnz);

      }
    }
    unit_assert(unit_incomplete);
    
    unit_func ("FluxData::set_neighbor_fluxes()");
  
    unit_assert(unit_incomplete);

    unit_func ("FluxData::delete_block_fluxes()");

    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
        flux_data.delete_block_fluxes (axis,face,index_field);
      }
    }
  
    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
        unit_assert(flux_data.block_fluxes(axis,face,index_field) == nullptr);
      }
    }
  
    unit_func ("FluxData::delete_neighbor_fluxes()");

    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
        flux_data.delete_neighbor_fluxes (axis,face,index_field);
      }
    }
  
    for (int axis=0; axis<3; axis++) {
      for (int face=0; face<2; face++) {
        unit_assert
          (flux_data.neighbor_fluxes(axis,face,index_field) == nullptr);
      }
    }

  }
  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

