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

//----------------------------------------------------------------------

EnzoBlock::EnzoBlock (CkMigrateMessage *m)
  : CBase_EnzoBlock (m)
{
  // replace Block's State with EnzoState
  state_ = std::make_shared<EnzoState>(0, 0.0, 0.0, false);

  TRACE("CkMigrateMessage");
  // EnzoSimulation[0] counts migrated Blocks
  proxy_enzo_simulation[0].p_method_balance_check();
}

//----------------------------------------------------------------------

EnzoBlock::EnzoBlock( process_type ip_source,  MsgType msg_type)
  : CBase_EnzoBlock (ip_source, msg_type)
{
  // replace Block's State with EnzoState
  state_ = std::make_shared<EnzoState>(0, 0.0, 0.0, false);
#ifdef TRACE_BLOCK

  CkPrintf ("%d %p TRACE_BLOCK %s EnzoBlock(ip) msg_type %d\n",
            CkMyPe(),(void *)this,name(thisIndex).c_str(),int(msg_type));
  fflush(stdout);
#endif

  if (msg_type == MsgType::msg_check) {
    set_msg_check(enzo::simulation()->get_msg_check(thisIndex));
  } else if (msg_type == MsgType::msg_refine) {
    set_msg_refine(enzo::simulation()->get_msg_refine(thisIndex));
  }
}

//----------------------------------------------------------------------

void EnzoBlock::set_msg_check(EnzoMsgCheck * msg)
{
  performance_start_(perf_block);

  restart_set_data_(msg);
  initialize();
  Block::initialize();
  performance_stop_(perf_block);
}

//----------------------------------------------------------------------

void EnzoBlock::set_msg_refine(MsgRefine * msg)
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
}

//======================================================================

void EnzoBlock::write(FILE * fp) throw ()
{
  const int in = cello::index_static();

  // Grid

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

  // problem

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

  // the following check was originally located elsewhere, but we moved it here
  // for lack of a better place to put it. We might be able to remove this in
  // the future.
  //
  // It's unclear whether we should be using field.num_permanent() or a
  // different count of fields. In any case, this is improvement over the count
  // that was used previously in this check
  if (field.num_permanent() > MAX_NUMBER_OF_BARYON_FIELDS) {
    ERROR2 ("EnzoBlock::initialize",
            "MAX_NUMBER_OF_BARYON_FIELDS = %d is too small for %d fields",
            MAX_NUMBER_OF_BARYON_FIELDS, field.num_permanent() );
  }

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  // query the ghost depth.
  //
  // Note: this is an improper way to do things... We only do it this way
  // because this is how it used to be done when EnzoBlock::ghost_depth was a
  // global variable
  //
  // In the future, we should be computing GridDimension, GridStartIndex and
  // GridEndIndex within the Method objects where that information is needed
  // and using the ghost_depth information relevant to the fields that are
  // being used. Specifically, this would look something like:
  //     int gx,gy,gz;
  //     field.ghost_depth(field.field_id("density"), &gx, &gy, &gz);
  //     if (cello::rank() < 1) gx = 0;
  //     if (cello::rank() < 2) gy = 0;
  //     if (cello::rank() < 3) gz = 0;

  const int rank = cello::rank();
  int gx = (rank < 1) ? 0 : cello::config()->field_ghost_depth[0];
  int gy = (rank < 2) ? 0 : cello::config()->field_ghost_depth[1];
  int gz = (rank < 3) ? 0 : cello::config()->field_ghost_depth[2];

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

bool EnzoBlock::spawn_child_blocks() throw()
{
  int level = index_.level();
  if (level >= 0) {
    
    if (level + 1 <= enzo::config()->refined_regions_lower.size()) {

      std::vector<int> lower = enzo::config()->refined_regions_lower.at(level);
      std::vector<int> upper = enzo::config()->refined_regions_upper.at(level);

      int ix, iy, iz, nx, ny, nz;
      index_global(&ix, &iy, &iz, &nx, &ny, &nz);
      if (lower.at(0) <= ix && ix < upper.at(0)) {
        if (lower.at(1) <= iy && iy < upper.at(1)) {
          if (lower.at(2) <= iz && iz < upper.at(2)) {
            return true;
          }
        }
      }
    }
  }

  return false;
}

void EnzoBlock::create_initial_child_blocks()
{
  bool spawn_children = spawn_child_blocks();
  if (spawn_children) {instantiate_children();}
}

void EnzoBlock::instantiate_children() throw()
{
  child_face_level_curr_.resize(cello::num_children()*27);
  int num_field_blocks = 1;

  int nx, ny, nz;
  data()->field().size(&nx, &ny, &nz);

  const int rank = cello::rank();
  ItChild it_child(rank);
  int ic3[3];
  // WARNING: code duplication in control_adapt.cpp adapt_refine_()
  while (it_child.next(ic3)) {
    Index index_child = index_.index_child(ic3);
    DataMsg * data_msg = NULL;

    std::vector<int> face_level;
    face_level.resize(27);
    const int o = 27*IC3(ic3);
    for (int i=0; i<face_level.size(); i++) {
      face_level[i] = child_face_level_curr_[o+i];
    }
    MsgRefine * msg = new MsgRefine
      (index_child,
       nx,ny,nz,
       num_field_blocks,
       adapt_step_,
       refresh_fine,
       face_level,
       &adapt_,state_.get());

    msg->set_data_msg(data_msg);
    cello::simulation()->p_refine_create_block (msg);

    children_.push_back(index_child);
  }

  adapt_.set_valid(false);
  is_leaf_ = false;
}
