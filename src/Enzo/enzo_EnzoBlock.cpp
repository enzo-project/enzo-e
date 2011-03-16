// $Id: enzo_EnzoBlock.cpp 2035 2011-02-28 23:47:31Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:02:02 PST 2011
/// @brief    Implementation of the EnzoBlock class

#include "cello.hpp"

#include "enzo.hpp"

//======================================================================

EnzoBlock::EnzoBlock(Patch * patch,
		     FieldDescr * field_descr,
		     int nx, int ny, int nz,
		     int num_field_blocks) throw()
  : Block(patch,field_descr,nx,ny,nz,num_field_blocks),
    CycleNumber(0),
    Time(0),
    OldTime(0),
    dt(0),
    SubgridFluxes(0)
{

  printf ("EnzoBlock()\n");
  int i,j;

  for (i=0; i<MAX_DIMENSION; i++) {
    AccelerationField[i] = 0;

    for (i=0; i<MAX_DIMENSION; i++) {
      AccelerationField[i] = 0;
      GridLeftEdge[i] = 0;
      GridDimension[i] = 0;
      GridStartIndex[i] = 0;
      GridEndIndex[i] = 0;
      CellWidth[i] = 0;
    }

    for (j=0; j<MAX_NUMBER_OF_BARYON_FIELDS; j++) {
      BaryonField[j] = 0;
      OldBaryonField[j] = 0;
    }

  }
  // CANNOT BE INITIALIZED HERE SINCE IT REQUIRES EXTENTS
  //  initialize();
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

EnzoBlock::~EnzoBlock() throw ()
{
}

//----------------------------------------------------------------------

void EnzoBlock::write(FILE * fp) throw ()
{
  fprintf (fp,"EnzoBlock: ComovingCoordinates %d\n",
	   ComovingCoordinates);
  fprintf (fp,"EnzoBlock: UseMinimumPressureSupport %d\n",
	   UseMinimumPressureSupport);
  fprintf (fp,"EnzoBlock: MinimumPressureSupportParameter %g\n",
	   MinimumPressureSupportParameter);
  fprintf (fp,"EnzoBlock: ComovingBoxSize %g\n",
	   ComovingBoxSize);
  fprintf (fp,"EnzoBlock: HubbleConstantNow %g\n",
	   HubbleConstantNow);
  fprintf (fp,"EnzoBlock: OmegaLambdaNow %g\n",
	   OmegaLambdaNow);
  fprintf (fp,"EnzoBlock: OmegaMatterNow %g\n",
	   OmegaMatterNow);
  fprintf (fp,"EnzoBlock: MaxExpansionRate %g\n",
	   MaxExpansionRate);

  // Chemistry

  fprintf (fp,"EnzoBlock: MultiSpecies %d\n",
	   MultiSpecies);

  // Gravity

  fprintf (fp,"EnzoBlock: GravityOn %d\n",
	   GravityOn);
  //  fprintf (fp,"EnzoBlock: *AccelerationField %g\n",
  //           *AccelerationField[MAX_DIMENSION)];

  // Physics

  fprintf (fp,"EnzoBlock: PressureFree %d\n",
	   PressureFree);
  fprintf (fp,"EnzoBlock: Gamma %g\n",
	   Gamma);
  fprintf (fp,"EnzoBlock: GravitationalConstant %g\n",
	   GravitationalConstant);

  // Problem-specific

  fprintf (fp,"EnzoBlock: ProblemType %d\n",
	   ProblemType);

  // Method PPM

  fprintf (fp,"EnzoBlock: PPMFlatteningParameter %d\n",
	   PPMFlatteningParameter);
  fprintf (fp,"EnzoBlock: PPMDiffusionParameter %d\n",
	   PPMDiffusionParameter);
  fprintf (fp,"EnzoBlock: PPMSteepeningParameter %d\n",
	   PPMSteepeningParameter);

  // Parallel

  fprintf (fp,"EnzoBlock: ProcessorNumber %d\n",
	   ProcessorNumber);

  // Numerics

  fprintf (fp,"EnzoBlock: DualEnergyFormalism %d\n",
	   DualEnergyFormalism);
  fprintf (fp,"EnzoBlock: DualEnergyFormalismEta1 %g\n",
	   DualEnergyFormalismEta1);
  fprintf (fp,"EnzoBlock: DualEnergyFormalismEta2 %g\n",
	   DualEnergyFormalismEta2);
  fprintf (fp,"EnzoBlock: pressure_floor %g\n",
	   pressure_floor);
  fprintf (fp,"EnzoBlock: density_density_floor %g\n",
	   density_floor);
  fprintf (fp,"EnzoBlock: number_density_floor %g\n",
	   number_density_floor);
  fprintf (fp,"EnzoBlock: temperature_floor %g\n",
	   temperature_floor);

  fprintf (fp,"EnzoBlock: CourantSafetyNumber %g\n",
	   CourantSafetyNumber);
  fprintf (fp,"EnzoBlock: InitialRedshift %g\n",
	   InitialRedshift);
  fprintf (fp,"EnzoBlock: InitialTimeInCodeUnits %g\n",
	   InitialTimeInCodeUnits);
  fprintf (fp,"EnzoBlock: Time %g\n",
	   Time);
  fprintf (fp,"EnzoBlock: OldTime %g\n",
	   OldTime);

  // Domain

  fprintf (fp,"EnzoBlock: DomainLeftEdge %g %g %g\n",
	   DomainLeftEdge [0],DomainLeftEdge [0],DomainLeftEdge [0]);
  fprintf (fp,"EnzoBlock: DomainRightEdge %g %g %g\n",
	   DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2]);

  // Fields

  if (field_density != -1) 
    fprintf (fp,"EnzoBlock: field_density %d\n", field_density);
  if (field_total_energy != -1) 
    fprintf (fp,"EnzoBlock: field_total_energy %d\n", field_total_energy);
  if (field_internal_energy != -1) 
    fprintf (fp,"EnzoBlock: field_internal_energy %d\n", field_internal_energy);
  if (field_velocity_x != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_x %d\n", field_velocity_x);
  if (field_velocity_y != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_y %d\n", field_velocity_y);
  if (field_velocity_z != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_z %d\n", field_velocity_z);
  if (field_color != -1) 
    fprintf (fp,"EnzoBlock: field_color %d\n", field_color);

  if (field_magnetic_x != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_x %d\n", field_magnetic_x);
  if (field_magnetic_y != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_y %d\n", field_magnetic_y);
  if (field_magnetic_z != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_z %d\n", field_magnetic_z);

  if (field_density_xp != -1) 
    fprintf (fp,"EnzoBlock: field_density_xp %d\n", field_density_xp);
  if (field_velocity_x_xp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_x_xp %d\n", field_velocity_x_xp);
  if (field_velocity_y_xp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_y_xp %d\n", field_velocity_y_xp);
  if (field_velocity_z_xp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_z_xp %d\n", field_velocity_z_xp);
  if (field_magnetic_x_xp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_x_xp %d\n", field_magnetic_x_xp);
  if (field_magnetic_y_xp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_y_xp %d\n", field_magnetic_y_xp);
  if (field_magnetic_z_xp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_z_xp %d\n", field_magnetic_z_xp);

  if (field_density_yp != -1) 
    fprintf (fp,"EnzoBlock: field_density_yp %d\n", field_density_yp);
  if (field_velocity_x_yp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_x_yp %d\n", field_velocity_x_yp);
  if (field_velocity_y_yp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_y_yp %d\n", field_velocity_y_yp);
  if (field_velocity_z_yp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_z_yp %d\n", field_velocity_z_yp);
  if (field_magnetic_x_yp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_x_yp %d\n", field_magnetic_x_yp);
  if (field_magnetic_y_yp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_y_yp %d\n", field_magnetic_y_yp);
  if (field_magnetic_z_yp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_z_yp %d\n", field_magnetic_z_yp);

  if (field_density_zp != -1) 
    fprintf (fp,"EnzoBlock: field_density_zp %d\n", field_density_zp);
  if (field_velocity_x_zp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_x_zp %d\n", field_velocity_x_zp);
  if (field_velocity_y_zp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_y_zp %d\n", field_velocity_y_zp);
  if (field_velocity_z_zp != -1) 
    fprintf (fp,"EnzoBlock: field_velocity_z_zp %d\n", field_velocity_z_zp);
  if (field_magnetic_x_zp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_x_zp %d\n", field_magnetic_x_zp);
  if (field_magnetic_y_zp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_y_zp %d\n", field_magnetic_y_zp);
  if (field_magnetic_z_zp != -1) 
    fprintf (fp,"EnzoBlock: field_magnetic_z_zp %d\n", field_magnetic_z_zp);

  // Grid

  fprintf (fp,"EnzoBlock: GridRank %d\n",    GridRank);
  fprintf (fp,"EnzoBlock: GridDimension %d %d %d\n",
	   GridDimension[0],GridDimension[1],GridDimension[2]);
  fprintf (fp,"EnzoBlock: GridStartIndex %d %d %d\n",
	   GridStartIndex[0],GridStartIndex[1],GridStartIndex[2]);
  fprintf (fp,"EnzoBlock: GridEndIndex %d %d %d\n",
	   GridEndIndex[0],GridEndIndex[1],GridEndIndex[2]);
  fprintf (fp,"EnzoBlock: GridLeftEdge %g %g %g\n",
	   GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);

  fprintf (fp,"EnzoBlock: CellWidth %g %g %g\n", 
	   CellWidth[0], CellWidth[1], CellWidth[2] );

  fprintf (fp,"EnzoBlock: ghost %d %d %d\n",
	   ghost_depth[0],ghost_depth[1],ghost_depth[2]);


  fprintf (fp,"EnzoBlock: NumberOfBaryonFields %d\n",
	   NumberOfBaryonFields);
  int i;
  for (i=0; i<NumberOfBaryonFields; i++) {
    fprintf (fp,"EnzoBlock: BaryonField[%d] %p\n",
	     i, BaryonField[i]);
    fprintf (fp,"EnzoBlock: OldBaryonField[%d] %p\n",
	     i, OldBaryonField[i]);
    fprintf (fp,"EnzoBlock: FieldType[%d] %d\n",
	     i, FieldType[i]);
  }

  fprintf (fp,"EnzoBlock: BoundaryRank %d\n", BoundaryRank);
  fprintf (fp,"EnzoBlock: BoundaryDimension %d %d %d\n",
	   BoundaryDimension[0],BoundaryDimension[1],BoundaryDimension[2]);

  // unknown, reflecting, outflow, inflow, periodic
  //  const char * bc_string[] = 
  //    {"unknown", "reflecting", "outflow", "inflow", "periodic"};

  for (i=0; i<NumberOfBaryonFields; i++) {

    fprintf (fp,"EnzoBlock: BoundaryFieldType[%d] %d\n", 
	     i, BoundaryFieldType[i]);

    fprintf (fp,"EnzoBlock: BoundaryType[%d] %p %p %p %p %p %p\n", i, 
	     BoundaryType[i][0][0],
	     BoundaryType[i][0][1],
	     BoundaryType[i][1][0],
	     BoundaryType[i][1][1],
	     BoundaryType[i][2][0],
	     BoundaryType[i][2][1]);

    fprintf (fp,"EnzoBlock: BoundaryValue[%d] %p %p %p %p %p %p\n",i,
	     BoundaryValue[i][0][0],
	     BoundaryValue[i][0][1],
	     BoundaryValue[i][1][0],
	     BoundaryValue[i][1][1],
	     BoundaryValue[i][2][0],
	     BoundaryValue[i][2][1]);  
  }

  // problem

  fprintf (fp,"EnzoBlock: CycleNumber %d\n",   CycleNumber);
  fprintf (fp,"EnzoBlock: dt %g\n", dt);

  // fluxes

  fprintf (fp,"EnzoBlock: SubgridFluxes %p\n", SubgridFluxes);
  

}

//----------------------------------------------------------------------
void EnzoBlock::initialize () throw()
{
  double xm,xp,ym,yp,zm,zp;

  lower(&xm,&ym,&zm);
  upper(&xp,&yp,&zp);

  GridLeftEdge[0]    = xm;
  GridLeftEdge[1]    = ym;
  GridLeftEdge[2]    = zm;

  // Grid dimensions

  int nx,ny,nz;
  field_block_[0] -> size (&nx,&ny,&nz);

  GridDimension[0]  = nx + 2*enzo::ghost_depth[0];
  GridDimension[1]  = ny + 2*enzo::ghost_depth[1];
  GridDimension[2]  = nz + 2*enzo::ghost_depth[2];
  GridStartIndex[0] = enzo::ghost_depth[0];
  GridStartIndex[1] = enzo::ghost_depth[1];
  GridStartIndex[2] = enzo::ghost_depth[2];
  GridEndIndex[0]   = enzo::ghost_depth[0] + nx - 1;
  GridEndIndex[1]   = enzo::ghost_depth[1] + ny - 1;
  GridEndIndex[2]   = enzo::ghost_depth[2] + nz - 1;

  // Initialize CellWidth

  double h3[3];
  field_block_[0]->cell_width(this,&h3[0],&h3[1],&h3[2]);

  for (int dim=0; dim<enzo::GridRank; dim++) {
    CellWidth[dim] = h3[dim];
  }

  // Initialize BaryonField[] pointers

  for (int field = 0; field < enzo::NumberOfBaryonFields; field++) {
    BaryonField[field] = (enzo_float *)field_block_[0]->field_values(field);
  }
}
