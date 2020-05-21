// See LICENSE_CELLO file for license and copyright information

/// @file     control_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with refreshing ghost zones
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define DEBUG_REFRESH

#ifdef DEBUG_REFRESH
#  define TRACE_REFRESH(msg,REFRESH)				\
  printf ("%d %s:%d %s TRACE_REFRESH %s type %d\n",CkMyPe(),	\
	  __FILE__,__LINE__,name().c_str(),msg,REFRESH->sync_type());	\
  fflush(stdout);
#else
#  define TRACE_REFRESH(msg,REFRESH) /* NOTHING */
#endif

//----------------------------------------------------------------------

void Block::refresh_begin_() 
{
  Refresh * refresh = this->refresh();
  TRACE_REFRESH("refresh_begin_()",refresh);

  check_leaf_();

  check_delete_();

  cello::simulation()->set_phase(phase_refresh);

  ERROR("refresh_begin_()",
        "refresh_begin_ called with NEW_REFRESH");
}

//----------------------------------------------------------------------

void Block::refresh_continue()
{

  // Refresh if Refresh object exists and have data
  Refresh * refresh = this->refresh();
  TRACE_REFRESH("refresh_continue_()",refresh);

  if ( refresh && refresh->is_active() ) {

    // count self
    int count = 1;

    // send Field face data
    if (refresh->any_fields()) {
      count += refresh_load_field_faces_ (refresh);
    }

    // send Particle face data
    if (refresh->any_particles()){
      count += refresh_load_particle_faces_ (refresh);
    }

    // wait for all messages to arrive (including this one)
    // before continuing to p_refresh_exit()

    ERROR("refresh_continue_()",
        "refresh_continue_ called with NEW_REFRESH");

  } else {

    refresh_exit_();
    
  }
  
}

//----------------------------------------------------------------------


void Block::p_refresh_store (MsgRefresh * msg)
{

  
  performance_start_(perf_refresh_store);

  msg->update(data());

  delete msg;

  Refresh * refresh = this->refresh();
  TRACE_REFRESH("p_refresh_store()",refresh);

  ERROR("p_refresh_store()",
        "p_refresh_store() called with NEW_REFRESH");
  
  performance_stop_(perf_refresh_store);
  performance_start_(perf_refresh_store_sync);
}


//----------------------------------------------------------------------

int Block::refresh_load_field_faces_ (Refresh *refresh)
{

  int count = 0;

  const int min_face_rank = refresh->min_face_rank();
  const int neighbor_type = refresh->neighbor_type();

  if (neighbor_type == neighbor_leaf ||
      neighbor_type == neighbor_tree) {

    // Loop over neighbor leaf Blocks (not necessarily same level)

    const int min_level = cello::config()->mesh_min_level;

    ItNeighbor it_neighbor =
      this->it_neighbor(min_face_rank,index_,
			neighbor_type,min_level,refresh->root_level());

    int if3[3];
    while (it_neighbor.next(if3)) {

      Index index_neighbor = it_neighbor.index();

      int ic3[3];
      it_neighbor.child(ic3);

      const int level = this->level();
      const int level_face = it_neighbor.face_level();

      const int refresh_type = 
	(level_face == level - 1) ? refresh_coarse :
	(level_face == level)     ? refresh_same :
	(level_face == level + 1) ? refresh_fine : refresh_unknown;

      refresh_load_field_face_ (refresh_type,index_neighbor,if3,ic3);
      ++count;
    }

  } else if (neighbor_type == neighbor_level) {

    // Loop over neighbor Blocks in same level (not necessarily leaves)

    ItFace it_face = this->it_face(min_face_rank,index_);

    int if3[3];
    while (it_face.next(if3)) {

      // count all faces if not a leaf, else don't count if face level
      // is less than this block's level
      
      if ( ! is_leaf() || face_level(if3) >= level()) {
	
	Index index_face = it_face.index();
	int ic3[3] = {0,0,0};
	refresh_load_field_face_ (refresh_same,index_face,if3,ic3);
	++count;
      }

    }
  }

  return count;
}

//----------------------------------------------------------------------

void Block::refresh_load_field_face_
( int refresh_type,
  Index index_neighbor,
  int if3[3],
  int ic3[3])

{
  //  TRACE_REFRESH("refresh_load_field_face()");

  // REFRESH FIELDS

  // ... coarse neighbor requires child index of self in parent

  if (refresh_type == refresh_coarse) {
    index_.child(index_.level(),ic3,ic3+1,ic3+2);
  }

  // ... copy field ghosts to array using FieldFace object
  bool lg3[3] = {false,false,false};

  Refresh * refresh = this->refresh();

  FieldFace * field_face = create_face
    (if3, ic3, lg3, refresh_type, refresh,false);
#ifdef DEBUG_FIELD_FACE  
  CkPrintf ("%d %s:%d DEBUG_FIELD_FACE creating %p\n",CkMyPe(),__FILE__,__LINE__,field_face);
#endif

  DataMsg * data_msg = new DataMsg;

  data_msg -> set_field_face (field_face,true);
  data_msg -> set_field_data (data()->field_data(),false);

  MsgRefresh * msg = new MsgRefresh;

  msg->set_data_msg (data_msg);

  thisProxy[index_neighbor].p_refresh_store (msg);

}


//----------------------------------------------------------------------

int Block::refresh_load_particle_faces_ (Refresh * refresh)
{

  //  TRACE_REFRESH("refresh_load_particle_faces()");
  
  const int rank = cello::rank();

  const int npa = (rank == 1) ? 4 : ((rank == 2) ? 4*4 : 4*4*4);

  ParticleData * particle_array[npa];
  ParticleData * particle_list [npa];
  Index * index_list = new Index[npa];
  
  for (int i=0; i<npa; i++) {
    particle_list[i]  = NULL;
    particle_array[i] = NULL;
  }

  // Sort particles that have left the Block into 4x4x4 array
  // corresponding to neighbors

  int nl = particle_load_faces_
    (npa,particle_list,particle_array, index_list, refresh);

  // Send particle data to neighbors

  particle_send_(nl,index_list,particle_list);

  delete [] index_list;

  return nl;
}

//----------------------------------------------------------------------

void Block::particle_send_
(int nl,Index index_list[], ParticleData * particle_list[])
{

  ParticleDescr * p_descr = cello::particle_descr();

  for (int il=0; il<nl; il++) {

    Index index           = index_list[il];
    ParticleData * p_data = particle_list[il];
    Particle particle_send (p_descr,p_data);
    
    if (p_data && p_data->num_particles(p_descr)>0) {

      DataMsg * data_msg = new DataMsg;
      data_msg ->set_particle_data(p_data,true);

      MsgRefresh * msg = new MsgRefresh;
      msg->set_data_msg (data_msg);

      thisProxy[index].p_refresh_store (msg);

    } else if (p_data) {
      
      MsgRefresh * msg = new MsgRefresh;

      thisProxy[index].p_refresh_store (msg);

      // assert ParticleData object exits but has no particles
      delete p_data;

    }

  }
}
