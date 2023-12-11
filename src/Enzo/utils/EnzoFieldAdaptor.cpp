// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFieldAdaptor.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-04-01
/// @brief    Implements the EnzoFieldAdaptor class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

// implement methods of BlockWrapper
namespace enzo_field_adaptor_detail {

  BlockWrapper::BlockWrapper(Block* block, int index_history)
    : block_(block),
      field_(block->data()->field()),
      index_history_(index_history)
  {
    ASSERT("BlockWrapper", "ghost zones must be allocated",
           field_.ghosts_allocated() == true);
  }

  //----------------------------------------------------------------------

  std::array<int,3> BlockWrapper::field_strides() const noexcept {
    int gx,gy;
    field_.ghost_depth (0,&gx,&gy,nullptr);

    int nx,ny;
    field_.size (&nx,&ny,nullptr);

    int ngx = nx + 2*gx;
    int ngy = ny + 2*gy;

    ASSERT("BlockWrapper::field_strides", "ghost zones must be allocated",
           field_.ghosts_allocated() == true); // sanity check

    return {ngx*ngy, ngx, 1};
  }

  //----------------------------------------------------------------------
  
  void BlockWrapper::grackle_field_grid_props(std::array<int,3>& grid_dimension,
                                              std::array<int,3>& grid_start,
                                              std::array<int,3>& grid_end)
    const noexcept
  {
    int gx,gy,gz;
    field_.ghost_depth (0,&gx,&gy,&gz);
    int nx,ny,nz;
    field_.size (&nx,&ny,&nz);
    int ngx = nx + 2*gx;
    int ngy = ny + 2*gy;
    int ngz = nz + 2*gz;

    grid_dimension = {ngx,     ngy,     ngz};

    // this includes the ghost zones
    grid_start   = {0,     0,     0};
    grid_end     = {ngx-1, ngy-1, ngz-1};
  }

  //----------------------------------------------------------------------

  double BlockWrapper::compute_time() const noexcept
  {
    if (index_history_ == 0) {
      return block_->time();
    } else {
      return field_.history_time(index_history_);
    }
  }

}

//----------------------------------------------------------------------

// implement methods of ArrayMapWrapper
namespace enzo_field_adaptor_detail {

  ArrayMapWrapper::ArrayMapWrapper(const EnzoEFltArrayMap& array_map)
    : array_map_(array_map)
  {
    ASSERT("ArrayMapWrapper", "array_map must hold at least 1 entry",
           array_map.size() >= 1);
  }

  //----------------------------------------------------------------------

  std::array<int,3> ArrayMapWrapper::field_strides() const noexcept{
    std::array<int,3> out;
    CelloView<const enzo_float,3> first_arr = array_map_[0];
    out = {first_arr.stride(0), first_arr.stride(1), first_arr.stride(2)};

    // TODO: track if strides are consistent within EnzoEFltArrayMap
    for (std::size_t i = 1; i < array_map_.size(); i++){
      CelloView<const enzo_float,3> cur = array_map_[i];
      if ( (out[0] != cur.stride(0)) || (out[1] != cur.stride(1)) ||
           (out[2] != cur.stride(2)) ){
        ERROR("ArrayMapWrapper::field_strides",
              "Not all of the contained arrays share strides");
      }
    }

    return out;
  }

  //----------------------------------------------------------------------

  void ArrayMapWrapper::grackle_field_grid_props
  (std::array<int,3>& grid_dimension,
   std::array<int,3>& grid_start,
   std::array<int,3>& grid_end) const noexcept
  {
    int nx = array_map_.array_shape(2);
    int ny = array_map_.array_shape(1);
    int nz = array_map_.array_shape(0);

    std::array<int,3> strides = this->field_strides();
    int stride_x = strides[2];
    // if we allow for selecting different indexing schemes in CelloView we'll
    // need to revise this method. In this scenario, strides[2] != 1
    ASSERT("ArrayMapWrapper::grackle_field_grid_props",
           "implmentation needs to be revisited", stride_x == 1);
    int stride_y = strides[1];
    int stride_z = strides[0];

    // stride_x == 1,  stride_y == ndx,  stride_z == ndx * ndy
    int ndx = stride_y;
    int ndy = stride_z / ndx;
    ASSERT("ArrayMapWrapper::grackle_field_grid_props",
           "stride_z % stride_y should be 0", stride_z % ndx == 0);
    int ndz = nz; // ndz is unimportant

    ASSERT("ArrayMapWrapper::grackle_field_grid_props", "sanity check",
           ( (ndx >= nx) & (ndy >= ny) & (ndz >= nz) ));

    grid_start = {0,0,0};
    grid_end = {nx - 1, ny - 1, nz - 1};
    grid_dimension = {ndx, ndy, ndz}; // if anything is wrong, it's this line

    if ((ndx != nx) | (ndy != ny)){
      ERROR("ArrayMapWrapper::try_get_ptr_grackle", "untested case");
    }
  }
}
