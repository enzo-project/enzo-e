// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoReconstructorNN.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Wed May 1 2019
/// @brief    [\ref Enzo] Implements the EnzoReconstructorNN class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoReconstructorNN::reconstruct_interface
(const EnzoEFltArrayMap &prim_map, EnzoEFltArrayMap &priml_map,
 EnzoEFltArrayMap &primr_map, const int dim, const EnzoEquationOfState *eos,
 const int stale_depth, const str_vec_t& passive_list)
{

  // determine components of i unit vector
  EnzoPermutedCoordinates coord(dim);
  int i_x, i_y, i_z;
  coord.i_unit_vector(i_z,i_y,i_x);

  auto fn = [i_x, i_y, i_z, coord, stale_depth,
             &prim_map, &priml_map, &primr_map](const std::string &key)
    {
      // define wc_offset(k,j,i) -> wc(k,j,i+1)
      const CelloArray<const enzo_float,3> wc = prim_map.get(key, stale_depth);
      const CelloArray<const enzo_float,3> wc_offset
        = coord.left_edge_offset(wc, 0, 0, 1);

      const EFlt3DArray wl = priml_map.get(key, stale_depth);
      const EFlt3DArray wr = primr_map.get(key, stale_depth);

      // the following could be optimized
      for (int iz=0; iz<wc.shape(0)-i_z; iz++) {
	for (int iy=0; iy<wc.shape(1)-i_y; iy++) {
	  for (int ix=0; ix<wc.shape(2)-i_x; ix++) {
	      wl(iz,iy,ix) = wc(iz,iy,ix);
	      wr(iz,iy,ix) = wc_offset(iz,iy,ix);
	  }
	}
      }
    };

  for (const std::string &key : active_key_names_){ fn(key); }
  for (const std::string &key : passive_list){ fn(key); }
}
