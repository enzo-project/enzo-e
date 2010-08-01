// $Id: user_MethodEnzoControl.cpp 1388 2010-04-20 23:57:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     user_MethodEnzoControl.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Create specific class for interfacing Cello code with User code
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of MethodEnzoControl user-dependent class member functions

#include "user.hpp"
#include "error.hpp"
#include "parameters.hpp"
#include "parallel.hpp"

#include "cello_hydro.h"

//----------------------------------------------------------------------

void MethodEnzoControl::initialize (DataDescr * data_descr) throw()
{
  // Extract field descriptor from data descriptor

  FieldDescr * field_descr = data_descr->field_descr();

  // Cosmology parameters

  Parameters * parameters = global_ -> parameters();

  parameters->set_current_group ("Physics");

  enzo_->ComovingCoordinates  = parameters->value_logical ("cosmology",false);
  Gamma                = parameters->value_scalar  ("gamma",5.0/3.0);

  //  if (ComovingCoordinates) {

    parameters->set_current_subgroup ("cosmology");

    InitialRedshift   = parameters->value_scalar ("initial_redshift",  20.0);
    HubbleConstantNow = parameters->value_scalar ("hubble_constant_now",0.701);
    OmegaLambdaNow    = parameters->value_scalar ("omega_lambda_now",   0.721);
    OmegaMatterNow    = parameters->value_scalar ("omega_matter_now",   0.279);
    MaxExpansionRate  = parameters->value_scalar ("max_expansion_rate", 0.01);
    ComovingBoxSize   = parameters->value_scalar ("comoving_box_size", 64.0);
    //  }

  int k = 0;

  field_density         = field_descr->insert_field("density");
  k = MAX(k,field_density);
  FieldType[field_density]      = Density;

  field_total_energy    = field_descr->insert_field("total_energy");
  k = MAX(k,field_total_energy);
  FieldType[field_total_energy] = TotalEnergy;

  if (DualEnergyFormalism) {
    field_internal_energy = field_descr->insert_field("internal_energy");
    k = MAX(k,field_internal_energy);
    FieldType[field_internal_energy] = InternalEnergy;
  }    

  parameters->set_current_group("Physics");
  GridRank = parameters->value_integer ("dimensions",0);

  if (GridRank >= 1) {
    field_velocity_x      = field_descr->insert_field("velocity_x");
    k = MAX(k,field_velocity_x);
    FieldType[field_velocity_x] = Velocity1;
  }

  if (GridRank >= 2) {
    field_velocity_y      = field_descr->insert_field("velocity_y");
    k = MAX(k,field_velocity_y);
    FieldType[field_velocity_y] = Velocity2;
  }

  if (GridRank >= 3) {
    field_velocity_z      = field_descr->insert_field("velocity_z");
    k = MAX(k,field_velocity_z);
    FieldType[field_velocity_z] = Velocity3;
  }

  field_color      = field_descr->insert_field("electron_density");
  k = MAX(k,field_color);
  FieldType[field_color] = ElectronDensity;

  // PPM parameters

  parameters->set_current_group ("Method","ppm");

  PressureFree = parameters->value_scalar("pressure_free",false);
  UseMinimumPressureSupport 
    =              parameters->value_logical("use_minimum_pressure_support",false);
  MinimumPressureSupportParameter 
    =              parameters->value_integer("minimum_pressure_support_parameter",100);
  PPMFlatteningParameter = parameters->value_logical ("flattening", false);
  PPMDiffusionParameter  = parameters->value_logical ("diffusion",  false);
  PPMSteepeningParameter = parameters->value_logical ("steepening", false);

  double floor_default = 1e-6;
  pressure_floor       = parameters->value_scalar("pressure_floor",      floor_default);
  density_floor        = parameters->value_scalar("density_floor",       floor_default);
  temperature_floor    = parameters->value_scalar("temperature_floor",   floor_default);
  number_density_floor = parameters->value_scalar("number_density_floor",floor_default);

  DualEnergyFormalism     = parameters->value_logical ("dual_energy",false);
  DualEnergyFormalismEta1 = parameters->value_scalar  ("dual_energy_eta_1",0.001);
  DualEnergyFormalismEta2 = parameters->value_scalar  ("dual_energy_eta_1",0.1);

  // Chemistry parameters

  MultiSpecies = 0;    // 0:0 1:6 2:9 3:12

  // Gravity parameters

  GravityOn                       = 0;    // Whether gravity is included
  GravitationalConstant           = 1.0;  // used only in SetMinimumSupport()
  AccelerationField[0]            = NULL;
  AccelerationField[1]            = NULL;
  AccelerationField[2]            = NULL;

  //Problem specific parameter

  ProblemType = 0;

  // Field parameters

  parameters->set_current_group ("Field");

  ghost_depth[0] = (GridRank >= 1) ? 
    parameters->list_value_integer(0,"ghosts",3) : 0;
  ghost_depth[1] = (GridRank >= 2) ? 
    parameters->list_value_integer(1,"ghosts",3) : 0;
  ghost_depth[2] = (GridRank >= 3) ? 
    parameters->list_value_integer(2,"ghosts",3) : 0;

  NumberOfBaryonFields = field_descr->field_count();

  ASSERT ("initialize",
	  "MAX_NUMBER_OF_BARYON_FIELDS is too small",
	  NumberOfBaryonFields <= MAX_NUMBER_OF_BARYON_FIELDS);

  Time                   = 0;

  // Domain parameters

  parameters->set_current_group ("Domain");
  
  DomainLeftEdge [0] = parameters->list_value_scalar(0,"extent",0.0);
  DomainRightEdge[0] = parameters->list_value_scalar(1,"extent",1.0);
  DomainLeftEdge [1] = parameters->list_value_scalar(2,"extent",0.0);
  DomainRightEdge[1] = parameters->list_value_scalar(3,"extent",1.0);
  DomainLeftEdge [2] = parameters->list_value_scalar(4,"extent",0.0);
  DomainRightEdge[2] = parameters->list_value_scalar(5,"extent",1.0);

  // Initial conditions

  parameters->set_current_group ("Initial");

  InitialTimeInCodeUnits = parameters->value_scalar ("time",0.0);
  Time = InitialTimeInCodeUnits;
  OldTime = Time;

  // Parallel parameters

  GroupProcessMpi parallel;

  ProcessorNumber = parallel.rank();

  parameters->set_current_group ("Field");
  
  CourantSafetyNumber = parameters->value_scalar ("courant",0.6);


}

//----------------------------------------------------------------------

void MethodEnzoControl::finalize (DataDescr * data_descr) throw()
{
}

//----------------------------------------------------------------------

void MethodEnzoControl::initialize_block ( DataBlock * data_block ) throw ()
{

  FieldBlock * field_block = data_block->field_block();
  
  double xm,xp,ym,yp,zm,zp;

  field_block->box_extent(&xm,&xp,&ym,&yp,&zm,&zp);

  GridLeftEdge[0]    = xm;
  GridLeftEdge[1]    = ym;
  GridLeftEdge[2]    = zm;

  // Grid dimensions

  int nx,ny,nz;
  field_block -> dimensions (&nx,&ny,&nz);

  GridDimension[0]  = nx + 2*ghost_depth[0];
  GridDimension[1]  = ny + 2*ghost_depth[1];
  GridDimension[2]  = nz + 2*ghost_depth[2];
  GridStartIndex[0] = ghost_depth[0];
  GridStartIndex[1] = ghost_depth[1];
  GridStartIndex[2] = ghost_depth[2];
  GridEndIndex[0]   = nx + ghost_depth[0] - 1;
  GridEndIndex[1]   = ny + ghost_depth[1] - 1;
  GridEndIndex[2]   = nz + ghost_depth[2] - 1;

  // Initialize CellWidth[].  Should be converted to constants hx,hy,hz

  double h3[3];
  field_block->cell_width(&h3[0],&h3[1],&h3[2]);

  for (int dim=0; dim<GridRank; dim++) {
    CellWidth[dim] = new ENZO_FLOAT[GridDimension[dim]];
    for (int i=0; i<GridDimension[dim]; i++) {
      CellWidth[dim][i] = h3[dim];
    }
  }

  // Initialize BaryonField[] pointers

  for (int field = 0; field < NumberOfBaryonFields; field++) {
    BaryonField[field] = (float *)field_block->field_values(field);
  }

 
//   // Boundary
//   /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
//      set it to one. */
 
//   /* 1) Compute Courant condition for baryons. */
 
//    // boundary
 
//    BoundaryRank = 2;
//    BoundaryDimension[0] = GridDimension[0];
//    BoundaryDimension[1] = GridDimension[1];

//    for (int field=0; field<NumberOfBaryonFields; field++) {
//      BoundaryFieldType[field] = FieldType[field];
//      for (int dim = 0; dim < 3; dim++) {
//        for (int face = 0; face < 2; face++) {
//  	int n1 = GridDimension[(dim+1)%3];
//  	int n2 = GridDimension[(dim+2)%3];
//  	int size = n1*n2;
//  	BoundaryType [field][dim][face] = new bc_type [size];
//  	BoundaryValue[field][dim][face] = NULL;
//  	for (int i2 = 0; i2<n2; i2++) {
//  	  for (int i1 = 0; i1<n1; i1++) {
//  	    int i = i1 + n1*i2;
//  	    BoundaryType[field][dim][face][i] = bc_reflecting;
//  	  }
//  	}
//        }
//      }
//    }

}


//----------------------------------------------------------------------

void MethodEnzoControl::finalize_block ( DataBlock * data_block ) throw ()
{
  ++ CycleNumber;
  for (int dim=0; dim<GridRank; dim++) {
    delete [] CellWidth[dim];
  }
}

//----------------------------------------------------------------------

void MethodEnzoControl::refresh_ghost(DataBlock * data_block, 
				      bool xm, bool xp, 
				      bool ym, bool yp, 
				      bool zm, bool zp) throw()
{

}

//======================================================================

