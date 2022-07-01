// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:02:02 PST 2011
/// @brief    Implementation of the EnzoBlock class

#include "cello.hpp"
#include "charm_simulation.hpp"

#include "enzo.hpp"

// #define TRACE_BLOCK

// #define DEBUG_ENZO_BLOCK

//======================================================================

int EnzoBlock::UseMinimumPressureSupport[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::MinimumPressureSupportParameter[CONFIG_NODE_SIZE];

// Chemistry

int EnzoBlock::MultiSpecies[CONFIG_NODE_SIZE];

// Physics

int EnzoBlock::PressureFree[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::GravitationalConstant[CONFIG_NODE_SIZE];

// Problem-specific

int EnzoBlock::ProblemType[CONFIG_NODE_SIZE];

// Method PPM

int EnzoBlock::PPMFlatteningParameter[CONFIG_NODE_SIZE];
int EnzoBlock::PPMDiffusionParameter[CONFIG_NODE_SIZE];
int EnzoBlock::PPMSteepeningParameter[CONFIG_NODE_SIZE];

// Numerics

enzo_float EnzoBlock::InitialRedshift[CONFIG_NODE_SIZE];
enzo_float EnzoBlock::InitialTimeInCodeUnits[CONFIG_NODE_SIZE];

// Domain

enzo_float EnzoBlock::DomainLeftEdge [3*CONFIG_NODE_SIZE];
enzo_float EnzoBlock::DomainRightEdge[3*CONFIG_NODE_SIZE];

// PPM

int EnzoBlock::GridRank[CONFIG_NODE_SIZE];

int EnzoBlock::ghost_depth[3*CONFIG_NODE_SIZE];

// Fields

int EnzoBlock::NumberOfBaryonFields[CONFIG_NODE_SIZE];

//----------------------------------------------------------------------

// STATIC
void EnzoBlock::initialize(const EnzoConfig * enzo_config)
{
#ifdef DEBUG_ENZO_BLOCK
  CkPrintf ("%d DEBUG_ENZO_BLOCK [static] EnzoBlock::initialize()\n",
            CkMyPe());
#endif
  int gx = enzo_config->field_ghost_depth[0];
  int gy = enzo_config->field_ghost_depth[1];
  int gz = enzo_config->field_ghost_depth[2];

  const int rank = enzo_config->mesh_root_rank;

  if (rank < 1) gx = 0;
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  double time  = enzo_config->initial_time;

  for (int in=0; in<CONFIG_NODE_SIZE; in++) {

    GridRank[in] = 0;
    NumberOfBaryonFields[in] = 0;

    int i;

    for (i=0; i<MAX_DIMENSION; i++) {
      DomainLeftEdge [in*3+i] = 0;
      DomainRightEdge[in*3+i] = 0;
      ghost_depth[in*3+i] = 0;
    }

    GridRank[in]            = enzo_config->mesh_root_rank;

    // Chemistry parameters

    MultiSpecies[in] = 0;    // 0:0 1:6 2:9 3:12

    // Gravity parameters

    GravitationalConstant[in]           = 1.0;  // used only in SetMinimumSupport()

    //Problem specific parameter

    ProblemType[in] = 0;

    // PPM parameters

    PressureFree[in]              = enzo_config->ppm_pressure_free;
    UseMinimumPressureSupport[in] = enzo_config->ppm_use_minimum_pressure_support;
    MinimumPressureSupportParameter[in] =
      enzo_config->ppm_minimum_pressure_support_parameter;

    PPMFlatteningParameter[in]    = enzo_config->ppm_flattening;
    PPMDiffusionParameter[in]     = enzo_config->ppm_diffusion;
    PPMSteepeningParameter[in]    = enzo_config->ppm_steepening;

    ghost_depth[in*3+0] = gx;
    ghost_depth[in*3+1] = gy;
    ghost_depth[in*3+2] = gz;

    NumberOfBaryonFields[in] = enzo_config->field_list.size();

    // Check NumberOfBaryonFields

    if (NumberOfBaryonFields[in] > MAX_NUMBER_OF_BARYON_FIELDS) {
      ERROR2 ("EnzoBlock::initialize",
	      "MAX_NUMBER_OF_BARYON_FIELDS = %d is too small for %d fields",
	      MAX_NUMBER_OF_BARYON_FIELDS,NumberOfBaryonFields[in] );
    }

    DomainLeftEdge [in*3+0] = enzo_config->domain_lower[0];
    DomainLeftEdge [in*3+1] = enzo_config->domain_lower[1];
    DomainLeftEdge [in*3+2] = enzo_config->domain_lower[2];

    DomainRightEdge[in*3+0] = enzo_config->domain_upper[0];
    DomainRightEdge[in*3+1] = enzo_config->domain_upper[1];
    DomainRightEdge[in*3+2] = enzo_config->domain_upper[2];

    InitialTimeInCodeUnits[in] = time;

  }

} // void initialize()

//----------------------------------------------------------------------

//======================================================================
#ifdef BYPASS_CHARM_MEM_LEAK
//======================================================================

EnzoBlock::EnzoBlock( process_type ip_source,  MsgType msg_type)
  : CBase_EnzoBlock (ip_source, msg_type),
    redshift(0.0)

{
#ifdef TRACE_BLOCK

  CkPrintf ("%d %p TRACE_BLOCK %s EnzoBlock(ip) msg_type %d\n",
            CkMyPe(),(void *)this,name(thisIndex).c_str(),int(msg_type));
#endif

  if (msg_type == MsgType::msg_check) {
    proxy_enzo_simulation[ip_source].p_get_msg_check(thisIndex);
  }
}

//----------------------------------------------------------------------

void EnzoBlock::p_set_msg_check(EnzoMsgCheck * msg)
{
  performance_start_(perf_block);

  restart_set_data_(msg);
  initialize();
  Block::initialize();
  performance_stop_(perf_block);
}

//----------------------------------------------------------------------

void EnzoBlock::p_set_msg_refine(MsgRefine * msg)
{
#ifdef TRACE_BLOCK
  CkPrintf ("%d %p :%d TRACE_BLOCK %s EnzoBlock p_set_msg_refine()\n",
            CkMyPe(),(void *)this,__LINE__,name(thisIndex).c_str());
  fflush(stdout);
#endif
  int io_reader = msg->restart_io_reader_;
  Block::p_set_msg_refine(msg);
  initialize();
  Block::initialize();
  // If refined block and restarting, notify file reader block is created
  if (io_reader >= 0) {
    proxy_io_enzo_reader[io_reader].p_block_created();
  }
}

//======================================================================
#else /* not BYPASS_CHARM_MEM_LEAK */
//======================================================================

EnzoBlock::EnzoBlock ( MsgRefine * msg )
  : CBase_EnzoBlock ( msg ),
    redshift(0.0)

{
#ifdef TRACE_BLOCK
  CkPrintf ("%d %p TRACE_BLOCK %s EnzoBlock(msg)\n",
            CkMyPe(),(void *)this,name(thisIndex).c_str());
#endif

  int io_reader = msg->restart_io_reader_;
  initialize();
  Block::initialize();
  // If refined block and restarting, notify file reader block is created
  if ((cello::config()->initial_restart) && (index_.level() > 0)) {
    proxy_io_enzo_reader[io_reader].p_block_created();
  }
  delete msg;
}

//----------------------------------------------------------------------

EnzoBlock::EnzoBlock ( EnzoMsgCheck * msg )
  : CBase_EnzoBlock (),
    redshift(0.0)
{
#ifdef TRACE_BLOCK
  CkPrintf ("%d %p TRACE_BLOCK %s EnzoBlock(msg)\n",
            CkMyPe(),(void *)this,name(thisIndex).c_str());
#endif

  init_refresh_();
  // init_refine_ (msg->index_,
  //       msg->nx_, msg->ny_, msg->nz_,
  //       msg->num_field_blocks_,
  //       msg->num_adapt_steps_,
  //       msg->cycle_, msg->time_,  msg->dt_,
  //       0, NULL, msg->refresh_type_,
  //       msg->num_face_level_, msg->face_level_,
  //       msg->adapt_parent_);

  // init_adapt_(msg->adapt_parent_);

  // apply_initial_(msg);

  initialize();
  Block::initialize();
#ifdef TRACE_BLOCK
  CkPrintf ("%d %p :%d TRACE_BLOCK %s EnzoBlock restart_set_data_\n",
            CkMyPe(),(void *)this,__LINE__,name(thisIndex).c_str());
  fflush(stdout);
#endif
  restart_set_data_(msg);
  delete msg;
}

//======================================================================
#endif
//======================================================================

EnzoBlock::~EnzoBlock()
{
#ifdef TRACE_BLOCK
  CkPrintf ("%d %p TRACE_BLOCK %s ~EnzoBlock(...)\n",
            CkMyPe(),(void *)this,name(thisIndex).c_str());
#endif
}

//----------------------------------------------------------------------

void EnzoBlock::pup(PUP::er &p)
{

  TRACEPUP;
  TRACE ("BEGIN EnzoBlock::pup()");

  CBase_EnzoBlock::pup(p);

  p | dt;

  const int in = cello::index_static();

  static bool warn1[CONFIG_NODE_SIZE] = {true};
  if (warn1[in]) {
    warn1[in] = false;
    WARNING("EnzoBlock::pup()", "skipping SubgridFluxes (not used)");
  }

  PUParray(p,GridLeftEdge,MAX_DIMENSION);
  PUParray(p,GridDimension,MAX_DIMENSION);
  PUParray(p,GridStartIndex,MAX_DIMENSION);
  PUParray(p,GridEndIndex,MAX_DIMENSION);
  PUParray(p,CellWidth,MAX_DIMENSION);

  p | redshift;
}

//======================================================================

void EnzoBlock::write(FILE * fp) throw ()
{
  const int in = cello::index_static();

  fprintf (fp,"EnzoBlock: UseMinimumPressureSupport %d\n",
	   UseMinimumPressureSupport[in]);
  fprintf (fp,"EnzoBlock: MinimumPressureSupportParameter %g\n",
	   MinimumPressureSupportParameter[in]);

  // Chemistry

  fprintf (fp,"EnzoBlock: MultiSpecies %d\n",
	   MultiSpecies[in]);

  // Physics

  fprintf (fp,"EnzoBlock: PressureFree %d\n",
	   PressureFree[in]);
  //fprintf (fp,"EnzoBlock: Gamma %g\n",
  //	   Gamma[in]);
  fprintf (fp,"EnzoBlock: GravitationalConstant %g\n",
	   GravitationalConstant[in]);

  // Problem-specific

  fprintf (fp,"EnzoBlock: ProblemType %d\n",
	   ProblemType[in]);

  // Method PPM

  fprintf (fp,"EnzoBlock: PPMFlatteningParameter %d\n",
	   PPMFlatteningParameter[in]);
  fprintf (fp,"EnzoBlock: PPMDiffusionParameter %d\n",
	   PPMDiffusionParameter[in]);
  fprintf (fp,"EnzoBlock: PPMSteepeningParameter %d\n",
	   PPMSteepeningParameter[in]);

  // Numerics

  /*
  fprintf (fp,"EnzoBlock: DualEnergyFormalism %d\n",
	   DualEnergyFormalism[in]);
  fprintf (fp,"EnzoBlock: DualEnergyFormalismEta1 %g\n",
	   DualEnergyFormalismEta1[in]);
  fprintf (fp,"EnzoBlock: DualEnergyFormalismEta2 %g\n",
	   DualEnergyFormalismEta2[in]);
  fprintf (fp,"EnzoBlock: pressure_floor %g\n",
	   pressure_floor[in]);
  fprintf (fp,"EnzoBlock: density_density_floor %g\n",
	   density_floor[in]);
  fprintf (fp,"EnzoBlock: number_density_floor %g\n",
	   number_density_floor[in]);
  fprintf (fp,"EnzoBlock: temperature_floor %g\n",
	   temperature_floor[in]);
  */

  fprintf (fp,"EnzoBlock: InitialRedshift %g\n",
	   InitialRedshift[in]);
  fprintf (fp,"EnzoBlock: InitialTimeInCodeUnits %g\n",
	   InitialTimeInCodeUnits[in]);

  // Domain

  fprintf (fp,"EnzoBlock: DomainLeftEdge %g %g %g\n",
	   DomainLeftEdge [in*3+0],
	   DomainLeftEdge [in*3+1],
	   DomainLeftEdge [in*3+2]);
  fprintf (fp,"EnzoBlock: DomainRightEdge %g %g %g\n",
	   DomainRightEdge[in*3+0],
	   DomainRightEdge[in*3+1],
	   DomainRightEdge[in*3+2]);

  // Fields

  // Grid

  fprintf (fp,"EnzoBlock: GridRank %d\n",    GridRank[in]);
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
	   ghost_depth[in*3+0],
	   ghost_depth[in*3+1],
	   ghost_depth[in*3+2]);


  fprintf (fp,"EnzoBlock: NumberOfBaryonFields %d\n",
	   NumberOfBaryonFields[in]);

  // problem

  fprintf (fp,"EnzoBlock: dt %g\n", dt);

}

//----------------------------------------------------------------------

void EnzoBlock::set_dt (double dt_param) throw ()
{
  Block::set_dt (dt_param);

  dt = dt_param;
}

//----------------------------------------------------------------------

void EnzoBlock::set_time (double time) throw ()
{
  Block::set_time (time);

  Simulation * simulation = cello::simulation();
  EnzoUnits * units = (EnzoUnits * )simulation->problem()->units();
  EnzoPhysicsCosmology * cosmology = units->cosmology();

  if (cosmology) {
    cosmology->set_current_time(time);
    redshift = cosmology->current_redshift();
  }
}

//----------------------------------------------------------------------

void EnzoBlock::initialize () throw()
{
  double xm,ym,zm;

  data()->lower(&xm,&ym,&zm);

  GridLeftEdge[0]  = xm;
  GridLeftEdge[1]  = ym;
  GridLeftEdge[2]  = zm;

  // Grid dimensions

  Field field = data()->field();

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int gx,gy,gz;

  const int in = cello::index_static();

  gx = EnzoBlock::ghost_depth[in*3+0];
  gy = EnzoBlock::ghost_depth[in*3+1];
  gz = EnzoBlock::ghost_depth[in*3+2];

  GridDimension[0]  = nx + 2*gx;
  GridDimension[1]  = ny + 2*gy;
  GridDimension[2]  = nz + 2*gz;

  TRACE("Initializing GridStartIndex");

  GridStartIndex[0] = gx;
  GridStartIndex[1] = gy;
  GridStartIndex[2] = gz;

  GridEndIndex[0] = gx + nx - 1;
  GridEndIndex[1] = gy + ny - 1;
  GridEndIndex[2] = gz + nz - 1;

  // Initialize CellWidth

  double xp,yp,zp;
  data()->upper(&xp,&yp,&zp);
  double hx,hy,hz;
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  CellWidth[0] = hx;
  CellWidth[1] = hy;
  CellWidth[2] = hz;

  TRACE ("Exit  EnzoBlock::initialize()\n");
}

