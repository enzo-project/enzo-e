// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#include "enzo.hpp"

#include "charm_enzo.hpp"

// #define TRACE_FACTORY

//----------------------------------------------------------------------

void EnzoFactory::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Factory::pup(p);
}

//----------------------------------------------------------------------

IoBlock * EnzoFactory::create_io_block () const throw()
{
  return new IoEnzoBlock;
}

//----------------------------------------------------------------------

CProxy_Block EnzoFactory::new_block_proxy
(
 DataMsg * data_msg,
 int nbx, int nby, int nbz) const throw()
{
  CProxy_EnzoBlock enzo_block_array;

  enzo_block_array = CProxy_EnzoBlock::ckNew();

  return enzo_block_array;
  
}

//----------------------------------------------------------------------
void EnzoFactory::create_block_array
(
 DataMsg * data_msg,
 CProxy_Block block_array,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks
 ) const throw()
{
#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY create_block_array %d %d %d\n",nbx,nby,nbz);
#endif  
  CProxy_EnzoBlock enzo_block_array = enzo::block_array();

  int count_adapt;

  int    cycle = enzo::simulation()->state()->cycle();
  double time  = enzo::simulation()->state()->time();
  double dt    = 0.0;
  State state (cycle,time,dt,false);

#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY %s:%d\n",__FILE__,__LINE__); fflush(stdout);
#endif

  int nax,nay,naz;
  cello::hierarchy()->root_blocks(&nax,&nay,&naz);

  std::vector<int> face_level;
  face_level.clear();

  for (int ix=0; ix<nbx; ix++) {
    for (int iy=0; iy<nby; iy++) {
      for (int iz=0; iz<nbz; iz++) {

        const int ip = ((long long) CkNumPes()) *
          (ix + nax*(iy + nay*iz)) / (nax*nay*naz);

        if (ip == CkMyPe()) {

          Index index(ix,iy,iz);

          MsgRefine * msg = new MsgRefine 
            (index,
             nx,ny,nz,
             num_field_blocks,
             count_adapt = 0,
             refresh_same,
             face_level,
             nullptr,
             &state);

          msg->set_data_msg(data_msg);

          // Use MappingArray initial mapping
          //  proxy_enzo_simulation[ip].p_refine_create_block (msg);
          enzo::simulation()->refine_create_block(msg);
        }
      }
    }
  }

  TRACE1("EnzoFactory::create_block_array = %p",&enzo_block_array);
}

//----------------------------------------------------------------------

void EnzoFactory::create_subblock_array
(
 DataMsg * data_msg,
 CProxy_Block block_array,
 int min_level,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 int num_field_blocks) const throw()
{
#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY create_subblock_array %d %d %d\n",nbx,nby,nbz);
#endif  
  TRACE8("EnzoFactory::create_subblock_array(min_level %d na(%d %d %d) n(%d %d %d) num_field_blocks %d",
	 min_level,nbx,nby,nbz,nx,ny,nz,num_field_blocks);

  if (min_level >= 0) {
    WARNING1("EnzoFactor::create_subblock_array",
	     "Trying to create subblock array with min_level %d >= 0",
	     min_level);
  }

  CProxy_EnzoBlock enzo_block_array = enzo::block_array();

  int nax,nay,naz;
  cello::hierarchy()->root_blocks(&nax,&nay,&naz);

  int    cycle = enzo::simulation()->state()->cycle();
  double time  = enzo::simulation()->state()->time();
  double dt    = 0.0;
  State state (cycle,time,dt,false);

  std::vector<int> face_level;
  face_level.clear();

  for (int level = -1; level >= min_level; level--) {

    if (nbx > 1) nbx = ceil(0.5*nbx);
    if (nby > 1) nby = ceil(0.5*nby);
    if (nbz > 1) nbz = ceil(0.5*nbz);

    int count_adapt;

#ifdef TRACE_FACTORY
    CkPrintf ("TRACE_FACTORY %s:%d\n",__FILE__,__LINE__); fflush(stdout);
#endif

    for (int ix=0; ix<nbx; ix++) {
      for (int iy=0; iy<nby; iy++) {
        for (int iz=0; iz<nbz; iz++) {

          const int ip = ((long long) CkNumPes()) *
            (ix + nax*(iy + nay*iz)) / (nax*nay*naz);

          if (ip == CkMyPe()) {

            int shift = -level;
 
            Index index(ix<<shift,iy<<shift,iz<<shift);

            index.set_level(level);

            TRACE3 ("inserting %d %d %d",ix,iy,iz);

            MsgRefine * msg = new MsgRefine 
              (index,
               nx,ny,nz,
               num_field_blocks,
               count_adapt=0,
               refresh_same,
               face_level,
               nullptr,
               &state);

            msg->set_data_msg(data_msg);

            enzo::simulation()->refine_create_block (msg);
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoFactory::create_block
(
 DataMsg * data_msg,
 CProxy_Block block_array,
 Index index,
 int nx, int ny, int nz,
 int num_field_blocks,
 int count_adapt,
 int narray, char * array, int refresh_type,
 const std::vector<int> & face_level,
 Adapt * adapt,
 State * state,
 Simulation * simulation,
 int io_reader,
 int ip) const throw()
{
#ifdef TRACE_FACTORY
  CkPrintf ("TRACE_FACTORY create_block %d %d %d\n");
  index.print("TRACE_FACTORY create_block",2);
#endif  

  TRACE3("EnzoFactory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("EnzoFactory::create_block() num_field_blocks %d  count_adapt %d",
	 num_field_blocks,count_adapt);

  CProxy_EnzoBlock enzo_block_array = (CProxy_EnzoBlock) block_array;

  const int rank = cello::rank();

  MsgRefine * msg = new MsgRefine 
    (index,
     nx,ny,nz,
     num_field_blocks,
     count_adapt,
     refresh_type,
     face_level,
     adapt,
     state,
     io_reader);

  msg->set_data_msg(data_msg);

  if (ip == -1) ip = CkMyPe();

  proxy_enzo_simulation[ip].p_refine_create_block (msg);
}

