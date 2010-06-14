// $Id: user_MethodEnzoPpm.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     user_MethodEnzoPpm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the MethodEnzoPpm class

#include "parallel.hpp"
#include "error.hpp"
#include "parameters.hpp"
#include "user_MethodEnzoPpm.hpp"

#include "ppm/cello_hydro.h"

void MethodEnzoPpm::initialize() throw()
{
  // Register method name

  method_name_ = "ppm";

  // Specify arguments

  add_argument_(argument_field, "density",        access_read_write);
  add_argument_(argument_field, "energy_total",   access_read_write);
  add_argument_(argument_field, "energy_internal",access_read_write);
  add_argument_(argument_field, "velocity_x",     access_read_write);
  add_argument_(argument_field, "velocity_y",     access_read_write);
  add_argument_(argument_field, "velocity_z",     access_read_write);

  // Initialize Enzo ppm hydrodynamics

  initialize_simulation_();

}

void MethodEnzoPpm::advance_block
(
 DataDescr * data_descr,
 DataBlock * data_block
 ) throw()
{
  INCOMPLETE_MESSAGE("MethodEnzoPpm::advance_block","");

  initialize_problem_(data_descr,data_block);

  // dt          = 
  // CycleNumber = 

  SolveHydroEquations (CycleNumber, dt);
}

void MethodEnzoPpm::refresh_face() throw()
{
  INCOMPLETE_MESSAGE("MethodEnzoPpm::refresh_face","");

}

void MethodEnzoPpm::initialize_simulation_ 
(
 DataDescr * data_descr
) throw()

{

  Parameters * p = Parameters::instance();

  p->set_current_group ("Method","ppm");


  // Cosmology parameters

  p->set_current_group ("Physics");

  GridRank             = p->value_integer ("dimensions",0);
  ComovingCoordinates  = p->value_logical ("cosmology",false);
  Gamma                = p->value_scalar  ("gamma",5.0/3.0);

  if (ComovingCoordinates) {

    p->set_current_subgroup ("cosmology");

    InitialRedshift   = p->value_scalar ("initial_redshift",  20.0);
    HubbleConstantNow = p->value_scalar ("hubble_constant_now",0.701);
    OmegaLambdaNow    = p->value_scalar ("omega_lambda_now",   0.721);
    OmegaMatterNow    = p->value_scalar ("omega_matter_now",   0.279);
    MaxExpansionRate  = p->value_scalar ("max_expansion_rate", 0.01);
    ComovingBoxSize   = p->value_scalar ("comoving_box_size", 64.0);
  }

  // Control parameters

  p->set_current_group ("Stopping");

  time_stop  = p->value_scalar ("time",0.0);
  cycle_stop = p->value_integer("cycle",0);

  // Field parameters

  p->set_current_group ("Field");
  
  CourantSafetyNumber    = p->value_scalar ("courant",0.6);

  int k = 0;

  WARNING_MESSAGE("MethodEnzoPpm::initialize_simulation_",
		  "Fixed fields");
  FieldType[field_density      = k++] = Density;
  FieldType[field_total_energy = k++] = TotalEnergy;
  FieldType[field_velocity_x   = k++] = Velocity1;
  FieldType[field_velocity_y   = k++] = Velocity2;
  FieldType[field_velocity_z   = k++] = Velocity3;

  // Domain parameters

  p->set_current_group ("Domain");
  
  DomainLeftEdge [0] = p->list_value_scalar(0,"extent",0.0);
  DomainRightEdge[0] = p->list_value_scalar(1,"extent",1.0);
  DomainLeftEdge [1] = p->list_value_scalar(2,"extent",0.0);
  DomainRightEdge[1] = p->list_value_scalar(3,"extent",1.0);
  DomainLeftEdge [2] = p->list_value_scalar(4,"extent",0.0);
  DomainRightEdge[2] = p->list_value_scalar(5,"extent",1.0);

  // Initial conditions

  p->set_current_group ("Initial");

  InitialTimeInCodeUnits = p->value_scalar ("time",0.0);

  // PPM parameters

  p->set_current_group ("Method","ppm");

  PressureFree              = p->value_scalar("pressure_free",false);
  UseMinimumPressureSupport 
    =              p->value_logical("use_minimum_pressure_support",false);
  MinimumPressureSupportParameter 
    =              p->value_integer("minimum_pressure_support_parameter",100);
  PPMFlatteningParameter = p->value_logical ("flattening", false);
  PPMDiffusionParameter  = p->value_logical ("diffusion",  false);
  PPMSteepeningParameter = p->value_logical ("steepening", false);

  double floor_default = 1e-20;
  pressure_floor       = p->value_scalar("pressure_floor",      floor_default);
  density_floor        = p->value_scalar("density_floor",       floor_default);
  temperature_floor    = p->value_scalar("temperature_floor",   floor_default);
  number_density_floor = p->value_scalar("number_density_floor",floor_default);

  DualEnergyFormalism     = p->value_logical ("dual_energy",false);
  DualEnergyFormalismEta1 = p->value_scalar  ("dual_energy_eta_1",0.001);
  DualEnergyFormalismEta2 = p->value_scalar  ("dual_energy_eta_1",0.1);

  // Chemistry parameters

  WARNING_MESSAGE("MethodEnzoPpm::initialize_simulation_",
		  "MultiSpecies parameters initialized to 0");
  
  MultiSpecies = 0;    // 0:0 1:6 2:9 3:12

  // Gravity parameters

  WARNING_MESSAGE("MethodEnzoPpm::initialize_simulation_",
		  "Gravity parameters initialized to 0");

  GravityOn                       = 0;    // Whether gravity is included
  GravitationalConstant           = 1.0;  // used only in SetMinimumSupport()
  AccelerationField[0]            = NULL;
  AccelerationField[1]            = NULL;
  AccelerationField[2]            = NULL;

  //Problem specific parameter

  WARNING_MESSAGE("MethodEnzoPpm::initialize_simulation_",
		  "ProblemType parameter set to 0");

  ProblemType = 0;

  // Parallel parameters

  ProcessorNumber = Parallel::instance()->process_rank();

  // Mesh parameters

  p->set_current_group ("Mesh");

  ghost_depth[0] = (GridRank >= 1) ? p->list_value_integer(0,"ghosts",3);
  ghost_depth[1] = (GridRank >= 2) ? p->list_value_integer(1,"ghosts",3);
  ghost_depth[2] = (GridRank >= 3) ? p->list_value_integer(2,"ghosts",3);

  // Boundary parameters

  INCOMPLETE_MESSAGE("MethodEnzoPpm::initialize_simulation_",
		     "Boundary parameters");

  BoundaryRank         = 0;
  BoundaryDimension[0] = 1;
  BoundaryDimension[1] = 1;
  BoundaryDimension[2] = 1;
 
  /* Clear BoundaryType and BoundaryValue pointers. */
 
  NumberOfBaryonFields = data_descr->field_count();

  for (int field = 0; field < num_fields; field++) {
    for (int dim = 0; dim < GridRank; dim++) {
      for (int i = 0; i < 2; i++) {
	BoundaryType[field][dim][i] = NULL;
	BoundaryValue[field][dim][i] = NULL;
      }
    }
  }
}

void MethodEnzoPpm::initialize_data_ 
(
 DataDescr * data_descr,
 DataBlock * data_block
 ) throw()
{
  FieldDescr * field_descr = data_descr->field_descr();
  FieldBlock * field_block = data_block->field_block();

  double xm,xp,ym,yp,zm,zp;
  field_block->box_extent(&xm,&ym,&zm,&xp,&yp,&zp);
  GridLeftEdge[0]    = xm;
  GridLeftEdge[1]    = ym;
  GridLeftEdge[2]    = zm;


}

void MethodEnzoPpm::initialize_problem_ 
(
 DataDescr * data_descr,
 DataBlock * data_block
) throw ()
{

  FieldBlock * field_block = data_block->field_block();
  
  // Get block dimensions

  int nx,ny,nz;
  field_block -> dimensions (&nx,&ny,&nz);

  CycleNumber = 0;
  dt          = 1;


  // Control

  Time                   = 0;
  OldTime                = 0;

  // Grid dimensions

  GridDimension[0]   = nx + 2*ghost_depth[0];
  GridDimension[1]   = ny + 2*ghost_depth[1];
  GridDimension[2]   = nz + 2*ghost_depth[2];
  GridStartIndex[0]  = ghost_depth[0];
  GridStartIndex[1]  = ghost_depth[1];
  GridStartIndex[2]  = ghost_depth[2];
  GridEndIndex[0]    = nx + ghost_depth[0] - 1;
  GridEndIndex[1]    = ny + ghost_depth[1] - 1;
  GridEndIndex[2]    = nz + ghost_depth[2] - 1;

  for (int dim=0; dim<GridRank; dim++) {
    CellWidth[dim] = new ENZO_FLOAT[GridDimension[dim]];
    float h = (DomainRightEdge[dim] - DomainLeftEdge[dim]) / 
      (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    for (int i=0; i<GridDimension[dim]; i++) {
      CellWidth[dim][i] = h;
    }
  }

  // Grid variables
  

  ASSERT ("initialize_implosion",
	  "MAX_NUMBER_OF_BARYON_FIELDS is too small",
	  NumberOfBaryonFields <= MAX_NUMBER_OF_BARYON_FIELDS);
  
  int nd = GridDimension[0] * GridDimension[1] * GridDimension[2];

  float * baryon_fields = new float [NumberOfBaryonFields * nd];
  for (int field = 0; field < NumberOfBaryonFields; field++) {
    BaryonField[field] = baryon_fields + field*nd;
  }

//   float * old_baryon_fields = new float [NumberOfBaryonFields * nd];
//   for (int field = 0; field < NumberOfBaryonFields; field++) {
//     OldBaryonField[field] = baryon_fields + field*nd;
//   }

  int ndx = GridDimension[0];
  //  int ndy = GridDimension[1];

  float hx = CellWidth[0][0];
  float hy = CellWidth[1][0];

  for (int iy = GridStartIndex[1]; iy<=GridEndIndex[1]; iy++) {

    float y = (iy - GridStartIndex[1] + 0.5)*hy;

    for (int ix = GridStartIndex[0]; ix<=GridEndIndex[0]; ix++) {

      float x = (ix - GridStartIndex[0] + 0.5)*hx;

      int i = ix + ndx * iy;

      // Initialize density

      if (x + y < cutoff) {
	BaryonField[ field_density ] [ i ] = density_in;
      } else {
	BaryonField[ field_density ] [ i ] = density_out;
      }

      // Initialize total energy

      if (x + y < cutoff) {
	BaryonField[ field_total_energy ][ i ] = 
	  pressure_in / ((Gamma - 1.0)*density_in);
      } else {
	BaryonField[ field_total_energy ][ i ] = 
	  pressure_out / ((Gamma - 1.0)*density_out);
      }

      // Initialize internal energy

      //      BaryonField[ field_internal_energy ][ i ] = 

      // Initialize velocity

      BaryonField[ field_velocity_x ][ i ] = velocity_x;
      BaryonField[ field_velocity_y ][ i ] = velocity_y;

      // Initialize color

      BaryonField[ field_color ][ i ] = 0.0;

    }
  }


  AccelerationField[0] = NULL;
  AccelerationField[1] = NULL;
  AccelerationField[2] = NULL;

  // Numerics

  DualEnergyFormalism             = 0;
  DualEnergyFormalismEta1         = 0.001;
  DualEnergyFormalismEta2         = 0.1;
  pressure_floor                  = 1e-6;
  number_density_floor            = 1e-6;
  density_floor                   = 1e-6;
  temperature_floor               = 1e-6;

  // boundary
 
  BoundaryRank = 2;
  BoundaryDimension[0] = GridDimension[0];
  BoundaryDimension[1] = GridDimension[1];

  for (int field=0; field<NumberOfBaryonFields; field++) {
    BoundaryFieldType[field] = FieldType[field];
    for (int dim = 0; dim < 3; dim++) {
      for (int face = 0; face < 2; face++) {
	int n1 = GridDimension[(dim+1)%3];
	int n2 = GridDimension[(dim+2)%3];
	int size = n1*n2;
	BoundaryType [field][dim][face] = new bc_type [size];
	BoundaryValue[field][dim][face] = NULL;
	for (int i2 = 0; i2<n2; i2++) {
	  for (int i1 = 0; i1<n1; i1++) {
	    int i = i1 + n1*i2;
	    BoundaryType[field][dim][face][i] = bc_reflecting;
	  }
	}
      }
    }
  }

