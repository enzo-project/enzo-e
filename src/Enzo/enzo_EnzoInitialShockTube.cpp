// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialShockTube.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2019-04-18
/// @brief    Implementation of EnzoInitialShockTube for initializing the RJ2a
///           shocktube problem mentioned by Stone et al 2008 and detailed in
///           Ryu & Jones (1995). This could probably be modified to include
///           additional shock tubes

// Note that there are a number of typos in various papers describing the
// initial conditions for the Ryu & Jones (95) 2a shock tube However, I think I
// have converged on the correct initial conditions from agreement between Dai
// & Woodward (93) (which lists the expected soln) and Stone et al. (08)
// Typos:
//   Ryu & Jones (95)
//     - states that the x-component of the magnetic field is of 2(4pi)^.5. It
//       should actually be 2/(4pi)^0.5
//     - They state that the left and right energies have values of 0.95 and 1.
//       Per Dai & Woodward (93) and Stone et al. (08) the left and right
//       pressures have values of 0.95 and 1.0.
//   Gardiner & Stone (08)
//     - states that both the left and right states have pressures of 0.95. Dai
//       & Woodward (93) and Stone et al. (08) indicate that the left state has 
//       a pressure of 0.95 and the right state has a pressure of 1.

#include "cello.hpp"
#include "enzo.hpp"
#include <math.h> // ceil

static std::vector<std::string> shock_tube_setups{"rj2a","sod"};



// right states for various shock tubes
static std::map<std::string, std::map<std::string, enzo_float>> shock_tube_l{
  {"rj2a", {// Ryu & Jones (95) from 2a
      {"density",    1.08},
      {"pressure",   0.95},
      {"velocity_0", 1.2},
      {"velocity_1", 0.01},
      {"velocity_2", 0.5},
      {"bfield_1",   1.0155412503859613},
      {"bfield_2",   0.5641895835477563}
  }},
  {"sod", {//Sod Shock Tube
      {"density",    1.0},
      {"pressure",   1.0},
      {"velocity_0", 0.},
      {"velocity_1", 0.},
      {"velocity_2", 0.},
      {"bfield_1",   0.},
      {"bfield_2",   0.}
  }}
};

// left states for different shock tubes
static std::map<std::string, std::map<std::string, enzo_float>> shock_tube_r{
  {"rj2a", {// Ryu & Jones (95) from 2a
      {"density", 1.},
      {"pressure", 1.0},
      {"velocity_0", 0.0},
      {"velocity_1", 0.0},
      {"velocity_2", 0.0},
      {"bfield_1", 1.1283791670955126},
      {"bfield_2", 0.5641895835477563}
  }},
  {"sod", {//Sod Shock Tube
      {"density",    0.125},
      {"pressure",   0.1},
      {"velocity_0", 0.},
      {"velocity_1", 0.},
      {"velocity_2", 0.},
      {"bfield_1",   0.},
      {"bfield_2",   0.}
  }}
};

// values bfield_0 for various shock tubes (aligned with the direction of
// propagation)
static std::map<std::string, enzo_float> shock_tube_bfield_0{
  {"rj2a", 0.5641895835477563},
  {"sod",  0.0}
};

//----------------------------------------------------------------------

std::string vector_to_string_(std::vector<std::string> &vec)
{
  std::string out = "";

  for(std::vector<std::string>::size_type i = 0; i != vec.size(); i++) {
    if (i != 0){
      out += ", ";
    }
    out += vec[i];
  }
  return out;
}

//----------------------------------------------------------------------

EnzoInitialShockTube::EnzoInitialShockTube(double gamma, int cycle, double time,
					   std::string setup_name,
					   std::string aligned_ax_name)
  : Initial(cycle, time), gamma_(gamma), setup_name_(setup_name), aligned_ax_(0)
{

  if (std::find(shock_tube_setups.begin(),
		shock_tube_setups.end(),
		setup_name_) == shock_tube_setups.end()){
    // the current name is not known
    std::string allowed_names = vector_to_string_(shock_tube_setups);
    // There is a character limit but we are probably fine (we are exiting
    // early anyways)
    ERROR2("EnzoInitialShockTube",
	   "Invalid setup_name specified (must be %s), not %s.",
	   allowed_names.c_str(), setup_name_.c_str());
  }


  ASSERT1("EnzoInitialShockTube",
	  "Invalid aligned_ax value specified (must be x, y, or z), not %s.",
	  aligned_ax_name.c_str(),
	  aligned_ax_name == "x" || aligned_ax_name == "y"
	  || aligned_ax_name == "z");
  if (aligned_ax_name == "x"){
    aligned_ax_ = 0;
  } else if (aligned_ax_name == "y"){
    aligned_ax_ = 1;
  } else {
    aligned_ax_ = 2;
  }
}
  
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
  std::map<std::string, enzo_float> l_vals = shock_tube_l[setup_name_];
  std::map<std::string, enzo_float> r_vals = shock_tube_r[setup_name_];
  enzo_float aligned_bfield_val = shock_tube_bfield_0[setup_name_];

  std::string velocities[3] = {"velocity_x", "velocity_y", "velocity_z"};
  std::string bfields[3] = {"bfieldi_x","bfieldi_y","bfieldi_z"};

  EnzoPermutedCoordinates coord(aligned_ax_);
  EnzoFieldArrayFactory array_factory(block);

  CSlice *l_slice = NULL;
  CSlice *r_slice = NULL;
  prep_aligned_slices_(block, &l_slice, &r_slice);

  for (int i = 0; i<2; i++){
    EFlt3DArray arr;
    CSlice *cur_slice = (i == 0) ? l_slice : r_slice;
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

    arr = array_factory.from_name(bfields[coord.j_axis()]);
    initializer_helper_(*cur_slice, cur_val_map->at("bfield_1"), arr);

    arr = array_factory.from_name(bfields[coord.k_axis()]);
    initializer_helper_(*cur_slice, cur_val_map->at("bfield_2"), arr);

    // compute the specific total energy
    arr = array_factory.from_name("total_energy");

    enzo_float etot, v2, b2;
    v2 = (cur_val_map->at("velocity_0") * cur_val_map->at("velocity_0") +
	  cur_val_map->at("velocity_1") * cur_val_map->at("velocity_1") +
	  cur_val_map->at("velocity_2") * cur_val_map->at("velocity_2"));
    b2 = (aligned_bfield_val * aligned_bfield_val +
	  cur_val_map->at("bfield_1") * cur_val_map->at("bfield_1") +
	  cur_val_map->at("bfield_2") * cur_val_map->at("bfield_2"));
    etot = ( (cur_val_map->at("pressure") /
	      ((gamma_ - 1.) * cur_val_map->at("density"))) +
	     0.5 * (v2 + b2 / cur_val_map->at("density")));

    initializer_helper_(*cur_slice, etot, arr);
  }

  EFlt3DArray align_b_arr = array_factory.from_name(bfields[coord.i_axis()]);
  align_b_arr.subarray() = aligned_bfield_val;

  delete l_slice;
  delete r_slice;

  // Compute the Cell-Centered B-fields
  EnzoInitialBCenter::initialize_bfield_center(block);
}

//----------------------------------------------------------------------

void EnzoInitialShockTube::prep_aligned_slices_(Block *block, CSlice **l_slice,
						CSlice **r_slice)
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
    *l_slice = new CSlice(0,shock_ind);
  }
  if (shock_ind < mi){
    *r_slice = new CSlice(shock_ind, mi);
  }
}

//----------------------------------------------------------------------

void EnzoInitialShockTube::initializer_helper_(CSlice &slice, enzo_float val,
					       EFlt3DArray &arr)
{
  CSlice xslice, yslice, zslice;
  
  xslice = (aligned_ax_ == 0) ? slice : CSlice(0, arr.shape(2));
  yslice = (aligned_ax_ == 1) ? slice : CSlice(0, arr.shape(1));
  zslice = (aligned_ax_ == 2) ? slice : CSlice(0, arr.shape(0));

  arr.subarray(zslice, yslice, xslice) = val;
}
