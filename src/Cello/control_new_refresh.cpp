// See LICENSE_CELLO file for license and copyright information

/// @file     control_new_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-05-23
/// @brief    Charm-related functions associated with refreshing ghost zones
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// #define DEBUG_NEW_REFRESH

#ifdef DEBUG_NEW_REFRESH
#   define TRACE_NEW_REFRESH(BLOCK,ID,MSG)			    \
  CkPrintf ("TRACE_NEW_REFRESH [%d] %s %s\n",ID,BLOCK->name().c_str(),MSG); \
  fflush(stdout);
#   define TRACE_STATE(BLOCK,MSG)					\
  CkPrintf ("TRACE_STATE %s %s state = %s\n",BLOCK->name().c_str(),MSG, \
	    (state==RefreshState::INACTIVE) ?				\
	    "inactive" : ((state == RefreshState::ACTIVE) ? "active" : "ready")); \
  fflush(stdout);
#else
#   define TRACE_NEW_REFRESH(BLOCK,ID,MSG) /* ... */
#   define TRACE_STATE(BLOCK,MSG) /* ... */
#endif

#define CHECK_ID(ID)				\
  ASSERT1 ("CHECK_ID",				\
	   "Invalid id %d",ID,			\
	   ID>=0);				\

#ifdef NEW_REFRESH

void Block::new_refresh_start (int id_refresh, int callback)
{
  CHECK_ID(id_refresh);
  TRACE_NEW_REFRESH(this,id_refresh,"new_refresh_start()");

  RefreshState & state = new_refresh_state_list_[id_refresh];
  Refresh & refresh    = new_refresh(id_refresh);

  // Send field and/or particle data associated with the given refresh
  // object to corresponding neighbors

#ifdef DEBUG_NEW_REFRESH  
    CkPrintf ("DEBUG_NEW_REFRESH %s is_active %d\n",name().c_str(),refresh.active());
#endif  
  
  if ( refresh.active() ) {

    ASSERT1 ("Block::new_refresh_start()",
	     "refresh[%d] state is not inactive",
	     id_refresh,
	     (state == RefreshState::INACTIVE));

    state = RefreshState::ACTIVE;
    TRACE_STATE(this,"new_refresh_start");

    int count = 0;

    // send Field face data
    
    if (refresh.any_fields()) {
      count += new_refresh_load_field_faces_ (refresh);
    }

    // send Particle face data
    if (refresh.any_particles()){
      count += new_refresh_load_particle_faces_(refresh);
    }

    Sync & sync = new_refresh_sync_list_[id_refresh];

    // Make sure sync counter is not active
    ASSERT4 ("Block::new_refresh_start()",
	     "refresh[%d] sync object %p is active (%d/%d)",
	     id_refresh, &sync, sync.value(), sync.stop(),
	     (sync.value() == 0 && sync.stop() == 0));
    
    // Initialize sync counter
    sync.set_stop(count);

#ifdef DEBUG_NEW_REFRESH
    CkPrintf ("DEBUG_NEW_REFRESH %s sync set %d/%d\n",
	      name().c_str(),sync.value(),sync.stop());
    fflush(stdout);
#endif    
    if (callback != 0) {
      new_refresh_wait(id_refresh,callback);
    }
  } else {
    if (callback != 0) 
      new_refresh_exit(refresh);
  }
}

//----------------------------------------------------------------------
void Block::new_refresh_wait (int id_refresh, int callback)
{
  CHECK_ID(id_refresh);
  TRACE_NEW_REFRESH(this,id_refresh,"new_refresh_wait()");

  Refresh & refresh = new_refresh(id_refresh);

  if (refresh.active()) {

    // make sure the callback parameter matches that in the refresh object

    ASSERT3("Block::new_refresh_wait()",
	   "Refresh[%d] mismatch between refresh callback %d and parameter callback %d",
	    id_refresh,callback, refresh.callback(),
	    (callback == refresh.callback()) );

    // make sure we aren't already in a "ready" state

    RefreshState & state = new_refresh_state_list_[id_refresh];
    
    ASSERT1("Block::new_refresh_wait()",
	   "Refresh[%d] not in 'active' state",
	    id_refresh,
	    (state == RefreshState::ACTIVE) );

    // tell refresh we're ready to start processing messages

    state = RefreshState::READY;
    TRACE_STATE(this,"new_refresh_wait");

    // process any existing messages in the refresh message list

    Sync & sync = new_refresh_sync_list_[id_refresh];

    for (auto id_msg=0;
	 id_msg<new_refresh_msg_list_[id_refresh].size();
	 id_msg++) {

      MsgRefresh * msg = new_refresh_msg_list_[id_refresh][id_msg];

      // unpack message data into Block data
      msg->update(data());
      
      delete msg;
      sync.next();
#ifdef DEBUG_NEW_REFRESH
    CkPrintf ("DEBUG_NEW_REFRESH %s sync next [active] %d/%d\n",
	      name().c_str(),sync.value(),sync.stop());
    fflush(stdout);
#endif    
    }

    // clear the message queue

    new_refresh_msg_list_[id_refresh].resize(0);

    // and check if we're finished

    new_refresh_check_done(id_refresh);
  }
}

//----------------------------------------------------------------------

void Block::new_refresh_check_done (int id_refresh)
{
  CHECK_ID(id_refresh);
  TRACE_NEW_REFRESH(this,id_refresh,"new_refresh_check_done()");

  Refresh & refresh    = new_refresh(id_refresh);
  RefreshState & state = new_refresh_state_list_[id_refresh];
  Sync & sync          = new_refresh_sync_list_[id_refresh];
  
  ASSERT1("Block::new_refresh_check_done()",
	  "Refresh[%d] must not be in inactive state",
	  id_refresh,
	  (state != RefreshState::INACTIVE) );

#ifdef DEBUG_NEW_REFRESH  
  CkPrintf ("DEBUG_NEW_REFRESH %s state==ready %d is_done %d %d/%d\n",
	    name().c_str(),(state == RefreshState::READY), sync.is_done(),
	    sync.value(),sync.stop());
#endif  
  
  TRACE_STATE(this,"new_refresh_check_done");
  
  if (sync.stop()==0 || state == RefreshState::READY && (sync.is_done())) {

    // Make sure incoming message queue is empty

    ASSERT2("Block::new_refresh_wait()",
	   "Refresh %d message list has size %d instead of 0",
	    id_refresh,new_refresh_msg_list_[id_refresh].size(),
	    (new_refresh_msg_list_[id_refresh].size() == 0));

    // reset sync counter
    sync.reset();
    sync.set_stop(0);
    
#ifdef DEBUG_NEW_REFRESH
    CkPrintf ("DEBUG_NEW_REFRESH %s sync reset %d/%d\n",
	      name().c_str(),sync.value(),sync.stop());
    fflush(stdout);
#endif    

    // reset refresh state to inactive

    state = RefreshState::INACTIVE;
    TRACE_STATE(this,"new_refresh_check_done");

    // Call callback

    new_refresh_exit(refresh);
  }
}

//----------------------------------------------------------------------

void Block::p_new_refresh_recv (MsgRefresh * msg)
{

  const int id_refresh = msg->id_refresh();
  CHECK_ID(id_refresh);
  TRACE_NEW_REFRESH(this,id_refresh,"p_new_refresh_recv()");
#ifdef DEBUG_NEW_REFRESH  
  CkPrintf ("DEBUG_NEW_REFRESH p_new_refresh_recv id %d\n",id_refresh);
#endif  

  RefreshState & state = new_refresh_state_list_[id_refresh];
  Sync & sync          = new_refresh_sync_list_[id_refresh];

  if (state == RefreshState::READY) {
    // unpack message data into Block data if ready
    msg->update(data());
      
    delete msg;
    sync.next();
#ifdef DEBUG_NEW_REFRESH
    CkPrintf ("DEBUG_NEW_REFRESH %s sync next [ready] %d/%d\n",
	      name().c_str(),sync.value(),sync.stop());
    fflush(stdout);
#endif    

    // check if it's the last message processed
    new_refresh_check_done(id_refresh);
  
  } else {

    // save message if not ready
    new_refresh_msg_list_[id_refresh].push_back(msg);

  }

}

//----------------------------------------------------------------------

void Block::new_refresh_exit (Refresh & refresh)
{
  TRACE_NEW_REFRESH(this,refresh.id(),"calling callback");
  CHECK_ID(refresh.id());
  update_boundary_();
  control_sync (refresh.callback(),
  		refresh.sync_type(),
  		refresh.sync_exit(),
  		refresh.min_face_rank(),
  		refresh.neighbor_type(),
  		refresh.root_level());
  // CkCallback
  //   (refresh.callback(),
  //    CkArrayIndexIndex(index_),thisProxy).send(NULL);
}

//----------------------------------------------------------------------

int Block::new_refresh_load_field_faces_ (Refresh & refresh)
{
  TRACE_NEW_REFRESH(this,refresh.id(),"new_refresh_load_field_faces_()");

  int count = 0;

  const int min_face_rank = refresh.min_face_rank();
  const int neighbor_type = refresh.neighbor_type();

  if (neighbor_type == neighbor_leaf ||
      neighbor_type == neighbor_tree) {

    // Loop over neighbor leaf Blocks (not necessarily same level)

    const int min_level = cello::config()->mesh_min_level;
    
    ItNeighbor it_neighbor =
      this->it_neighbor(min_face_rank,index_,
			neighbor_type,min_level,refresh.root_level());

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

      new_refresh_load_field_face_ (refresh,refresh_type,index_neighbor,if3,ic3);
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
	new_refresh_load_field_face_ (refresh,refresh_same,index_face,if3,ic3);
	++count;
      }

    }
  }

  return count;
}

//----------------------------------------------------------------------

void Block::new_refresh_load_field_face_
( Refresh & refresh,
  int refresh_type,
  Index index_neighbor,
  int if3[3],
  int ic3[3])

{
  TRACE_NEW_REFRESH(this,refresh.id(),"new_refresh_load_field_face_()");
  // ... coarse neighbor requires child index of self in parent

  if (refresh_type == refresh_coarse) {
    index_.child(index_.level(),ic3,ic3+1,ic3+2);
  }

  // ... copy field ghosts to array using FieldFace object

  MsgRefresh * msg_refresh = new MsgRefresh;

  DataMsg * data_msg = new DataMsg;

  bool lg3[3] = {false,false,false};

  FieldFace * field_face = create_face
    (if3, ic3, lg3, refresh_type, &refresh,false);

  data_msg -> set_field_face (field_face,true);
  data_msg -> set_field_data (data()->field_data(),false);

  const int id_refresh = refresh.id();
  CHECK_ID(id_refresh);
  
  ASSERT1 ("Block::new_refresh_load_field_face_()",
	  "id_refresh %d of refresh object is out of range",
	   id_refresh,
	   (0 <= id_refresh));
  msg_refresh->set_new_refresh_id (id_refresh);
  msg_refresh->set_data_msg (data_msg);

  thisProxy[index_neighbor].p_new_refresh_recv (msg_refresh);

}


//----------------------------------------------------------------------

int Block::new_refresh_load_particle_faces_ (Refresh & refresh)
{
  TRACE_NEW_REFRESH(this,refresh.id(),"new_refresh_load_particle_faces_()");
  const int rank = cello::rank();

  const int npa3[3] = { 4, 4*4, 4*4*4 };
  const int npa = npa3[rank-1];

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
    (npa,particle_list,particle_array, index_list, &refresh);

  // Send particle data to neighbors

  new_particle_send_(refresh,nl,index_list,particle_list);

  delete [] index_list;

  return nl;
}

//----------------------------------------------------------------------

void Block::new_particle_send_
(Refresh & refresh, int nl,Index index_list[], ParticleData * particle_list[])
{
  TRACE_NEW_REFRESH(this,-1,"new_particle_send_()");

  ParticleDescr * p_descr = cello::particle_descr();

  for (int il=0; il<nl; il++) {

    Index index           = index_list[il];
    ParticleData * p_data = particle_list[il];
    Particle particle_send (p_descr,p_data);
    
    const int id_refresh = refresh.id();
    CHECK_ID(id_refresh);

    ASSERT1 ("Block::new_refresh_load_field_face_()",
	     "id_refresh %d of refresh object is out of range",
	     id_refresh,
	     (0 <= id_refresh));

  if (p_data && p_data->num_particles(p_descr)>0) {

      DataMsg * data_msg = new DataMsg;
      data_msg ->set_particle_data(p_data,true);

      MsgRefresh * msg_refresh = new MsgRefresh;
      msg_refresh->set_data_msg (data_msg);
      msg_refresh->set_new_refresh_id (id_refresh);

      thisProxy[index].p_new_refresh_recv (msg_refresh);

    } else if (p_data) {
      
      MsgRefresh * msg_refresh = new MsgRefresh;

      msg_refresh->set_data_msg (nullptr);
      msg_refresh->set_new_refresh_id (id_refresh);

      thisProxy[index].p_new_refresh_recv (msg_refresh);

      // assert ParticleData object exits but has no particles
      delete p_data;

    }

  }
}
#endif // NEW_REFRESH
