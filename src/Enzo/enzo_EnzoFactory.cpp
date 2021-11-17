// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#include "enzo.hpp"

#include "charm_enzo.hpp"

// #define DEBUG_NEW_ADAPT

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

  //  CProxy_MappingTree array_map  = CProxy_MappingTree::ckNew(nbx,nby,nbz);
  CProxy_MappingArray array_map  = CProxy_MappingArray::ckNew(nbx,nby,nbz);


  CkArrayOptions opts;
  opts.setMap(array_map);

  enzo_block_array = CProxy_EnzoBlock::ckNew(opts);

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

  CProxy_EnzoBlock enzo_block_array = enzo::block_array();

  int count_adapt;

  int    cycle = 0;
  double time  = 0.0;
  double dt    = 0.0;
  int num_face_level = 0;
  int * face_level = 0;

  for (int ix=0; ix<nbx; ix++) {
    for (int iy=0; iy<nby; iy++) {
      for (int iz=0; iz<nbz; iz++) {

	Index index(ix,iy,iz);

	MsgRefine * msg = new MsgRefine 
	  (index,
	   nx,ny,nz,
	   num_field_blocks,
	   count_adapt = 0,
	   cycle,time,dt,
	   refresh_same,
	   num_face_level, face_level);

	msg->set_data_msg(data_msg);

	enzo::simulation()->set_msg_refine (index,msg);
	enzo_block_array[index].insert (process_type(CkMyPe()));

	// --------------------------------------------------

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
  TRACE8("EnzoFactory::create_subblock_array(min_level %d na(%d %d %d) n(%d %d %d) num_field_blocks %d",
	 min_level,nbx,nby,nbz,nx,ny,nz,num_field_blocks);

  if (min_level >= 0) {
    WARNING1("EnzoFactor::create_subblock_array",
	     "Trying to create subblock array with min_level %d >= 0",
	     min_level);
  }

  CProxy_EnzoBlock enzo_block_array = enzo::block_array();

  for (int level = -1; level >= min_level; level--) {

    if (nbx > 1) nbx = ceil(0.5*nbx);
    if (nby > 1) nby = ceil(0.5*nby);
    if (nbz > 1) nbz = ceil(0.5*nbz);

    // CProxy_MappingTree array_map  = CProxy_MappingTree::ckNew(nbx,nby,nbz);
    CProxy_MappingArray array_map  = CProxy_MappingArray::ckNew(nbx,nby,nbz);

    CkArrayOptions opts;
    opts.setMap(array_map);

    int count_adapt;

    int    cycle = 0;
    double time  = 0.0;
    double dt    = 0.0;
    int num_face_level = 0;
    int * face_level = 0;

    for (int ix=0; ix<nbx; ix++) {
      for (int iy=0; iy<nby; iy++) {
	for (int iz=0; iz<nbz; iz++) {

	  int shift = -level;
	  
	  Index index(ix<<shift,iy<<shift,iz<<shift);

	  index.set_level(level);

	  TRACE3 ("inserting %d %d %d",ix,iy,iz);

	  MsgRefine * msg = new MsgRefine 
	    (index,
	     nx,ny,nz,
	     num_field_blocks,
	     count_adapt=0,
	     cycle,time,dt,
	     refresh_same,
	     num_face_level, face_level);

	  msg->set_data_msg(data_msg);

	  enzo::simulation()->set_msg_refine (index,msg);
	  enzo_block_array[index].insert (process_type(CkMyPe()));

	  // --------------------------------------------------

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
 int cycle, double time, double dt,
 int narray, char * array, int refresh_type,
 int num_face_level,
 int * face_level,
 Adapt * adapt,
 Simulation * simulation
 ) const throw()
{

  TRACE3("EnzoFactory::create_block(%d %d %d)",nx,ny,nz);
  TRACE2("EnzoFactory::create_block() num_field_blocks %d  count_adapt %d",
	 num_field_blocks,count_adapt);

  CProxy_EnzoBlock enzo_block_array = (CProxy_EnzoBlock) block_array;

  const int rank = cello::rank();

  int iym=(rank >= 2) ? 0 : 1;
  int iyp=(rank >= 2) ? 3 : 2;
  int izm=(rank >= 3) ? 0 : 1;
  int izp=(rank >= 3) ? 3 : 2;
#ifdef NEW_ADAPT
  CkPrintf ("adapt %p\n",adapt); fflush(stdout);

  char buffer[80];
  int v3[3];
  index.values(v3);
  sprintf (buffer,"Adapt [ %8X %8X %8X ]\n",v3[0],v3[1],v3[2]);
  adapt->print(buffer);

  int * face_level_new = new int[27];
  std::fill_n(face_level_new,27,INDEX_UNDEFINED_LEVEL);

  for (int i=0; i<adapt->num_neighbors(); i++) {
    Index index_neighbor = adapt->index(i);
    int level = index_neighbor.level();
    int im3[3]={0,0,0},ip3[3]={1,1,1};
    index.categorize (index_neighbor,rank,im3,ip3);
    for (int iz=im3[2]; iz<ip3[2]; iz++) {
      for (int iy=im3[1]; iy<ip3[1]; iy++) {
        for (int ix=im3[0]; ix<ip3[0]; ix++) {
          const int i = ix + 3*(iy + 3*iz);
          face_level_new[i] = level;
        }
      }
    }
  }

  CkPrintf ("TRACE %d\n",__LINE__); fflush(stdout);

  // int level = index.level();
  // int cx,cy,cz;
  // index.child(level,&cx,&cy,&cz);
  // for (int iz=izm; iz<3; iz++) {
  //   int ipz=(iz+1+cz) / 2;
  //   for (int iy=iym; iy<3; iy++) {
  //     int ipy=(iy+1+cy) / 2;
  //     for (int ix=0; ix<3; ix++) {
  //       int ipx=(ix+1+cx) / 2;
  //       int i=ix+3*(iy+3*iz);
  //       int ip=ipx+3*(ipy+3*ipz);
  //       face_level[i] = face_level_parent[ip];
  //     }
  //   }
  // }
  // // Increment faces internal to block (including self) to be in fine level 
  // iyp=(rank >= 2) ? 2 : 1;
  // izp=(rank >= 3) ? 2 : 1;
  // for (int iz=0; iz<izp; iz++) {
  //   for (int iy=0; iy<iyp; iy++) {
  //     for (int ix=0; ix<2; ix++) {
  //       int i = (ix+1-cx) + 3*( (iy+1-cy) + 3*(iz+1-cz));
  //       ++ face_level[i];
  //     }
  //   }
  // }
#endif
#ifdef DEBUG_NEW_ADAPT  
  iym=(rank >= 2) ? 0 : 1;
  iyp=(rank >= 2) ? 3 : 2;
  izm=(rank >= 3) ? 0 : 1;
  izp=(rank >= 3) ? 3 : 2;
  int nb3[3] = {2,2,1};
  index.print("face level",2,2,nb3,true);
  for (int iz=izm; iz<izp; iz++) {
    for (int iy=iym; iy<iyp; iy++) {
      for (int ix=0; ix<3; ix++) {
        int i=ix+3*(iy+3*iz);
        CkPrintf ("%1d",face_level[i]);
      }
    }
  }
  CkPrintf ("\n");
#ifdef NEW_ADAPT
  // CkPrintf ("Child %d %d %d\n",cx,cy,cz);
  // index.print("face level parent",2,2,nb3,true);
  // for (int iz=izm; iz<izp; iz++) {
  //   for (int iy=iym; iy<iyp; iy++) {
  //     for (int ix=0; ix<3; ix++) {
  //       int i=ix+3*(iy+3*iz);
  //       CkPrintf ("%1d",face_level_parent[i]);
  //     }
  //   }
  // }
  // CkPrintf ("\n");
#endif  /* NEW_ADAPT */
  
#endif /* DEBUG_NEW!_ADAPT */
  MsgRefine * msg = new MsgRefine 
    (index,
     nx,ny,nz,
     num_field_blocks,
     count_adapt,
     cycle,time,dt,
     refresh_type,
     num_face_level, face_level);

  msg->set_data_msg(data_msg);

  enzo::simulation()->set_msg_refine (index,msg);

  enzo_block_array[index].insert ( process_type(CkMyPe()) );
}

