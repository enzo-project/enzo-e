// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShockTube.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-04-18
/// @brief    Implementation of EnzoInitialShockTube for initializing the RJ2a
///           shocktube problem mentioned by Stone et al 2008 and detailed in
///           Ryu & Jones (1995). This could probably be modified to include
///           additional shock tubes

#include "cello.hpp"
#include "enzo.hpp"
#include <math.h> // ceil

//----------------------------------------------------------------------

void EnzoInitialShockTube::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  p | aligned_ax_;
}

//----------------------------------------------------------------------

void EnzoInitialShockTube::enforce_block 
( Block * block, const Hierarchy  * hierarchy ) throw()

{
  std::map<std::string, enzo_float> l_vals = {
    {"density", 1.08}, {"pressure", 0.95},
    {"velocity_0", 1.2}, {"velocity_1", 0.01}, {"velocity_2", 0.5},
    {"bfield_1", 1.0155412503859613}, {"bfield_2", 0.5641895835477563}
  };

  std::map<std::string, enzo_float> r_vals = {
    {"density", 1.}, {"pressure", 0.95},
    {"velocity_0", 0.0}, {"velocity_1", 0.0}, {"velocity_2", 0.0},
    {"bfield_1", 1.1283791670955126}, {"bfield_2", 0.5641895835477563}
  };
  enzo_float aligned_bfield_val = 0.5641895835477563;

  std::string velocities[3] = {"velocity_x", "velocity_y", "velocity_z"};
  std::string bfields[3] = {"bfieldi_x","bfieldi_y","bfieldi_z"};

  EnzoPermutedCoordinates coord(aligned_ax_);
  EnzoFieldArrayFactory array_factory(block);

  ESlice *l_slice = NULL;
  ESlice *r_slice = NULL;
  prep_aligned_slices_(block, &l_slice, &r_slice);

  for (int i = 0; i<2; i++){
    EFlt3DArray arr;
    ESlice *cur_slice = (i == 0) ? l_slice : r_slice;
    std::map<std::string, enzo_float> *cur_val_map;
    cur_val_map = (i == 0) ? &l_vals : &r_vals;
    if (cur_slice == NULL){
      continue;
    }
    arr = array_factory.from_name("density");
    initializer_helper_(*cur_slice, cur_val_map->at("density"), arr);

    arr = array_factory.from_name(velocities[coord.i_axis()]);
    initializer_helper_(*cur_slice, cur_val_map->at("velocity_0"),arr);

    arr = array_factory.from_name(velocities[coord.j_axis()]);
    initializer_helper_(*cur_slice, cur_val_map->at("velocity_1"), arr);

    arr = array_factory.from_name(velocities[coord.k_axis()]);
    initializer_helper_(*cur_slice, cur_val_map->at("velocity_2"), arr);

    //arr = array_factory.from_name("pressure");
    //initializer_helper_(*cur_slice, cur_val_map->at("pressure"), arr);

    arr = array_factory.from_name(bfields[coord.j_axis()]);
    initializer_helper_(*cur_slice, cur_val_map->at("bfield_1"), arr);

    arr = array_factory.from_name(bfields[coord.k_axis()]);
    initializer_helper_(*cur_slice, cur_val_map->at("bfield_2"), arr);
  }
  EFlt3DArray align_b_arr = array_factory.from_name(bfields[coord.i_axis()]);
  align_b_arr.subarray() = aligned_bfield_val;

  EFlt3DArray pressure_array = array_factory.from_name("pressure");
  pressure_array.subarray() = 0.95;
  

  delete l_slice;
  delete r_slice;

  // Compute the Cell-Centered B-fields
  EnzoInitialBCenter::initialize_bfield_center(block);
}

//----------------------------------------------------------------------

void EnzoInitialShockTube::prep_aligned_slices_(Block *block, ESlice **l_slice,
						ESlice **r_slice)
{
  Field field  = block->data()->field();
  
  // total size of the block (including ghost zones)
  int mx,my,mz;
  field.dimensions (field.field_id("density"),&mx,&my,&mz);

  // ghost depth
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);
  
  // cell widths
  double xmb, ymb, zmb, xpb, ypb, zpb, hx, hy, hz;
  block->data()->lower(&xmb,&ymb,&zmb); // lower left edge of active region
  block->data()->upper(&xpb,&ypb,&zpb); // upper right edge of active region
  field.cell_width(xmb,xpb,&hx,
		   ymb,ypb,&hy,
		   zmb,zpb,&hz);

  // let axis i be the aligned axis
  int mi, gi;
  double di, left_edge;
  if (aligned_ax_ == 0){
    mi = mx;   gi = gx;   di = hx;   left_edge = xmb;
  } else if (aligned_ax_ == 1){
    mi = my;   gi = gy;   di = hy;   left_edge = ymb;
  } else {
    mi = mz;   gi = gz;   di = hz;   left_edge = zmb;
  }

  // the cell-centered position corresponding to index ind is:
  //   left_edge + di * (0.5+(double)(ind-gi))
  // The first cell with pos>=0
  int shock_ind = (int)ceil((0.5-left_edge)/di-0.5 + (double)gi);
  
  if (shock_ind > 0){
    *l_slice = new ESlice(0,shock_ind);
  }
  if (shock_ind < mi){
    *r_slice = new ESlice(shock_ind, mi);
  }
}

//----------------------------------------------------------------------

void EnzoInitialShockTube::initializer_helper_(ESlice &slice, enzo_float val,
					       EFlt3DArray &arr)
{
  ESlice xslice, yslice, zslice;
  
  xslice = (aligned_ax_ == 0) ? slice : ESlice(0, arr.shape(2));
  yslice = (aligned_ax_ == 1) ? slice : ESlice(0, arr.shape(1));
  zslice = (aligned_ax_ == 2) ? slice : ESlice(0, arr.shape(0));

  arr.subarray(zslice, yslice, xslice) = val;
}
