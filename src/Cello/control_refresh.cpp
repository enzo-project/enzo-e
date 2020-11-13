// See LICENSE_CELLO file for license and copyright information

/// @file     control_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-05-23
/// @brief    Charm-related functions associated with refreshing ghost zones
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#define CHECK_ID(ID) ASSERT1 ("CHECK_ID","Invalid id %d",ID,(ID>=0));

// #define DEBUG_ENZO_PROLONG
// #define TRACE_COMM

//----------------------------------------------------------------------

#ifdef TRACE_COMM

#   define TRACE_COMM_SEND(TYPE,BLOCK,NAME_RECV,ID)                     \
  {                                                                     \
    CkPrintf ("TRACE_COMM %s send %s >> %s send %d\n",   \
              TYPE,BLOCK->name().c_str(),NAME_RECV.c_str(),ID);         \
  }
#   define TRACE_COMM_RECV(BLOCK,ID,MSG_REFRESH)                        \
  {                                                                     \
    CkPrintf ("TRACE_COMM %s recv %s << %s recv %d\n",   \
              MSG_REFRESH->type_name().c_str(),                         \
              BLOCK->name().c_str(),                                    \
              MSG_REFRESH->block_name().c_str(), ID);                   \
  }

#   define TRACE_COMM_XPCT(TYPE,NAME_RECV,ID,NAME_SEND)                 \
  {                                                                     \
    CkPrintf ("TRACE_COMM %s xpct %s << %s xpct %d\n",   \
              TYPE,                                                     \
              NAME_RECV.c_str(),                                        \
              NAME_SEND.c_str(), ID);                                   \
  }

#else

#   define TRACE_COMM_SEND(MSG,BLOCK,NAME_RECV,ID) /* ... */
#   define TRACE_COMM_RECV(BLOCK,ID,MSG_REFRESH) /* ... */
#   define TRACE_COMM_XPCT(TYPE,BLOCK,ID,NAME_RECV) /* ... */

#endif

//----------------------------------------------------------------------

void Block::refresh_start (int id_refresh, int callback)
{
  CHECK_ID(id_refresh);
  Refresh * refresh = cello::refresh(id_refresh);
  Sync * sync = sync_(id_refresh);

  // Send field and/or particle data associated with the given refresh
  // object to corresponding neighbors

  if ( refresh->is_active() ) {

    ASSERT1 ("Block::refresh_start()",
	     "refresh[%d] state is not inactive",
	     id_refresh,
	     (sync->state() == RefreshState::INACTIVE));

    sync->set_state(RefreshState::ACTIVE);

    // send Field face data

    int count_field=0;
    if (refresh->any_fields()) {
      count_field = refresh_load_field_faces_ (*refresh);
    }

    // send Particle face data
    int count_particle=0;
    if (refresh->any_particles()){
      count_particle = refresh_load_particle_faces_(*refresh);
    }

    // send Flux face data
    int count_flux=0;
    if (refresh->any_fluxes()){
      count_flux = refresh_load_flux_faces_(*refresh);
    }

    const int count = count_field + count_particle + count_flux;
    
    // Make sure sync counter is not active
    ASSERT4 ("Block::refresh_start()",
	     "refresh[%d] sync object %p is active (%d/%d)",
	     id_refresh, sync, sync->value(), sync->stop(),
	     (sync->value() == 0 && sync->stop() == 0));

    // Initialize sync counter
    sync->set_stop(count);
#ifdef TRACE_COMM    
    CkPrintf ("TRACE_COMM %s %d count %d\n",name().c_str(),id_refresh,count);
#endif    

    refresh_wait(id_refresh,callback);

  } else {

    refresh_exit(*refresh);

  }
}

//----------------------------------------------------------------------
void Block::refresh_wait (int id_refresh, int callback)
{
  CHECK_ID(id_refresh);

  Refresh * refresh = cello::refresh(id_refresh);
  Sync * sync = sync_(id_refresh);

  ASSERT1("Block::refresh_wait()",
          "Wait called with inactive Refresh[%d]",
          id_refresh,
          (refresh->is_active()));

  // make sure the callback parameter matches that in the refresh object

  ASSERT3("Block::refresh_wait()",
          "Refresh[%d] mismatch between refresh %d and parameter %d callbacks",
          id_refresh,callback, refresh->callback(),
          (callback == refresh->callback()) );

  // make sure we aren't already in a "ready" state

  ASSERT1("Block::refresh_wait()",
          "Refresh[%d] not in 'active' state",
          id_refresh,
          (sync->state() == RefreshState::ACTIVE) );

  // tell refresh we're ready to start processing messages

  sync->set_state(RefreshState::READY);

  // process any existing messages in the refresh message list

  for (auto id_msg=0;
       id_msg<refresh_msg_list_[id_refresh].size();
       id_msg++) {

    MsgRefresh * msg = refresh_msg_list_[id_refresh][id_msg];

    // unpack message data into Block data
    msg->update(data());
      
    TRACE_COMM_RECV(this,id_refresh,msg);
    delete msg;
    sync->advance();
  }
    
  // clear the message queue

  refresh_msg_list_[id_refresh].resize(0);

  // and check if we're finished

  refresh_check_done(id_refresh);
}

//----------------------------------------------------------------------

void Block::refresh_check_done (int id_refresh)
{
  CHECK_ID(id_refresh);

  Refresh * refresh = cello::refresh(id_refresh);
  Sync * sync = sync_(id_refresh);

  ASSERT1("Block::refresh_check_done()",
	  "Refresh[%d] must not be in inactive state",
	  id_refresh,
	  (sync->state() != RefreshState::INACTIVE) );

  if ( (sync->stop()==0) ||
      (sync->is_done() && (sync->state() == RefreshState::READY))) {

    // Make sure incoming message queue is empty

    ASSERT2("Block::refresh_wait()",
	   "Refresh %d message list has size %lu instead of 0",
	    id_refresh,refresh_msg_list_[id_refresh].size(),
	    (refresh_msg_list_[id_refresh].size() == 0));

    // reset sync counter
    sync->reset();
    sync->set_stop(0);
    
    // reset refresh state to inactive

    sync->set_state(RefreshState::INACTIVE);

    // Call callback

    refresh_exit(*refresh);
  }
}

//----------------------------------------------------------------------

void Block::p_refresh_recv (MsgRefresh * msg)
{

  const int id_refresh = msg->id_refresh();
  CHECK_ID(id_refresh);

  Sync * sync = sync_(id_refresh);

  if (sync->state() == RefreshState::READY) {

    // unpack message data into Block data if ready
    msg->update(data());
      
    delete msg;

    sync->advance();
    TRACE_COMM_RECV(this,id_refresh,msg);

    // check if it's the last message processed
    refresh_check_done(id_refresh);
  
  } else {

    // save message if not ready
    refresh_msg_list_[id_refresh].push_back(msg);

  }

}

//----------------------------------------------------------------------

void Block::refresh_exit (Refresh & refresh)
{
  CHECK_ID(refresh.id());
  update_boundary_();
  control_sync (refresh.callback(),
  		refresh.sync_type(),
  		refresh.sync_exit(),
  		refresh.min_face_rank(),
  		refresh.neighbor_type(),
  		refresh.root_level());
}

//----------------------------------------------------------------------

int Block::refresh_load_field_faces_ (Refresh & refresh)
{
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

      // handle padded interpolation special case if needed

      count += refresh_load_extra_face_
        (refresh,index_neighbor, level,level_face,if3,ic3);

      refresh_load_field_face_ (refresh,refresh_type,index_neighbor,if3,ic3);
      ++count;
      TRACE_COMM_XPCT("field",name(),refresh.id(),name(index_neighbor));
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
	refresh_load_field_face_ (refresh,refresh_same,index_face,if3,ic3);
	++count;
        TRACE_COMM_XPCT("field",name(),refresh.id(),name(index_face));

      }

    }
  }

  return count;
}

//----------------------------------------------------------------------

void Block::refresh_load_field_face_
( Refresh & refresh,
  int refresh_type,
  Index index_neighbor,
  int if3[3],
  int ic3[3])

{
  // ... coarse neighbor requires child index of self in parent

  if (refresh_type == refresh_coarse) {
    index_.child(index_.level(),ic3,ic3+1,ic3+2);
  }

  // ... copy field ghosts to array using FieldFace object

  MsgRefresh * msg_refresh = new MsgRefresh;
#ifdef TRACE_MSG_REFRESH  
  msg_refresh->set_block_type(name(),"field");
#endif
  bool lg3[3] = {false,false,false};

  FieldFace * field_face = create_face
    (if3, ic3, lg3, refresh_type, &refresh,false);

  DataMsg * data_msg = new DataMsg;

  data_msg -> set_field_face (field_face,true);
  data_msg -> set_field_data (data()->field_data(),false);

  const int id_refresh = refresh.id();
  CHECK_ID(id_refresh);

  msg_refresh->set_refresh_id (id_refresh);
  msg_refresh->set_data_msg (data_msg);

  TRACE_COMM_SEND("field",this,name(index_neighbor),id_refresh);
  thisProxy[index_neighbor].p_refresh_recv (msg_refresh);

}

//----------------------------------------------------------------------

int Block::refresh_load_extra_face_
(Refresh refresh,
 Index index_neighbor,
 int level, int level_face,
 int if3[3],
 int ic3[3])
{
  int count = 0;

  const int padding = cello::problem()->prolong()->padding();

  if ((padding > 0) && (level != level_face)) {

    // Create box_face
    
    const int rank = cello::rank();
    int n3[3];
    data()->field().size(n3,n3+1,n3+2);
    const int g = refresh.ghost_depth();
    int g3[3] = {g3[0] = (rank >= 1) ? g : 0,
                 g3[1] = (rank >= 2) ? g : 0,
                 g3[2] = (rank >= 3) ? g : 0};

    // Boxes used Bs -> br, Be -> br, Bs -> be
    // (s send, r receive, e extra)
    Box box_sr (rank,n3,g3);
    Box box_se (rank,n3,g3);
    Box box_er (rank,n3,g3);
    
    // Create iterator over extra blocks
    
    ItNeighbor it_extra =
      this->it_neighbor(refresh.min_face_rank(),index_,
                        refresh.neighbor_type(),
                        cello::config()->mesh_min_level,
                        refresh.root_level());
    
      // ... determine intersection region

    const bool l_send = (level < level_face);
    const bool l_recv = (level > level_face);

    int jf3[3] = { l_send ? if3[0] : -if3[0],
                   l_send ? if3[1] : -if3[1],
                   l_send ? if3[2] : -if3[2] };
    
    box_sr.set_block(+1,jf3,ic3);
    box_sr.set_padding(padding);

    box_sr.compute_region();
#ifdef DEBUG_ENZO_PROLONG    
    if (l_send) box_sr.print("sr region send");
    if (l_recv) box_sr.print("sr region recv");
#endif    

    if (l_send) {

      // SENDER LOOP OVER EXTRA BLOCKS
      
      int ef3[3];
      while (it_extra.next(ef3)) {

        const Index index_extra = it_extra.index();
        const int   level_extra = it_extra.face_level();
        
        int ec3[3] = {0,0,0};
        if (level_extra > level) {
          index_extra.child(level_extra,ec3,ec3+1,ec3+2);
        }
#ifdef DEBUG_ENZO_PROLONG        
        CkPrintf ("DEBUG_EXTRA %s send loop: if3 %d %d %d ic3 %d %d %d\n",
                  name().c_str(),
                  if3[0],if3[1],if3[2],
                  ic3[0],ic3[1],ic3[2]);
        CkPrintf ("DEBUG_EXTRA %s send loop: ef3 %d %d %d ic3 %d %d %d\n",
                  name().c_str(),
                  ef3[0],ef3[1],ef3[2],
                  ec3[0],ec3[1],ec3[2]);
        CkPrintf ("DEBUG_EXTRA %s send loop: Ls Lr Le %d %d %d\n",
                  name().c_str(),
                  level,level_face,level_extra);
#endif        
        
        // ... skip extra block if it's the same as the neighbor
        
        bool l_valid_level = (std::abs(level_extra - level) <= 1);

        if (index_extra != index_neighbor && l_valid_level) {

#ifdef DEBUG_ENZO_PROLONG
          CkPrintf ("DEBUG_EXTRA %s send uniq\n",name().c_str());
#endif          
          // ... determine overlap of extra block with intersection region
          const int level_send = level;
          box_sr.set_block ((level_extra-level_send), ef3,ec3);
          box_sr.compute_block_start();
#ifdef DEBUG_ENZO_PROLONG
          box_sr.print("sr overlap send");
#endif          
        
          int im3[3],ip3[3];
          bool overlap = box_sr.get_limits (im3,ip3,Box::BlockType::extra);
          
#ifdef DEBUG_ENZO_PROLONG
          CkPrintf ("DEBUG_PROLONG send overlap %d if3 %d %d ef3 %d %d im3 %d %d ip3 %d %d\n",
                    overlap?1:0,
                    if3[0],if3[1],
                    ef3[0],ef3[1],
                    im3[0],im3[1],
                    ip3[0],ip3[1]);
#endif          
          if (overlap) {

#ifdef DEBUG_ENZO_PROLONG
            CkPrintf ("DEBUG_EXTRA %s send overlap\n",name().c_str());
#endif            
            if (level_extra == level) {
          
#ifdef DEBUG_ENZO_PROLONG
              CkPrintf ("DEBUG_EXTRA %s send extra\n",name().c_str());
#endif              
              // this block sends; extra block is coarse
          
              // handle contribution of this block Bs to Be -> br

              int if3_er[3] = { if3[0]-ef3[0], if3[1]-ef3[1], if3[2]-ef3[2] };
                                    

              // Box Bs | Be -> br
              box_er.set_block(+1,if3_er,ic3);
              box_er.set_padding(padding);
              box_er.compute_region();

              int itm3[3],itp3[3];
              box_er.get_limits (itm3,itp3,Box::BlockType::send);

              int if3_es[3] = {-ef3[0], -ef3[1], -ef3[2] };
              
              box_er.set_block(0,if3_es,ic3); // ic3 ignored
              box_er.compute_block_start();

              int ifm3[3],ifp3[3];
              int iam3[3],iap3[3];
              box_er.get_limits (ifm3,ifp3,Box::BlockType::extra);
              box_er.get_limits (iam3,iap3,Box::BlockType::array);

              int ma3[3] =
                { itp3[0]-itm3[0], itp3[1]-itm3[1], itp3[2]-itm3[2]};

              refresh_extra_send_
                (refresh, index_neighbor, if3_er,
                 ma3, iam3,iap3, ifm3,ifp3, data()->field());

              ASSERT9 ("Block::refresh_load_extra_face_",
                       "Array limits %d %d %d - %d %d %d not within array size %d %d %d\n",
                       iam3[0],iam3[1],iam3[2],
                       iap3[0],iap3[1],iap3[2],
                       ma3[0], ma3[1], ma3[2],
                       (0 <= iam3[0] && iap3[0] <= ma3[0]) &&
                       (0 <= iam3[1] && iap3[1] <= ma3[1]) &&
                       (0 <= iam3[2] && iap3[2] <= ma3[2]));
                       
#ifdef DEBUG_ENZO_PROLONG
              CkPrintf ("DEBUG_PROLONG send array size  %d %d %d\n",
                        ma3[0],ma3[1],ma3[2]);
              CkPrintf ("DEBUG_PROLONG send limits ifm3 %d %d %d ifp3 %d %d %d\n",
                        ifm3[0],ifm3[1],ifm3[2],ifp3[0],ifp3[1],ifp3[2]);
              CkPrintf ("DEBUG_PROLONG send limits iam3 %d %d %d iap3 %d %d %d\n\n",
                        iam3[0],iam3[1],iam3[2],iap3[0],iap3[1],iap3[2]);
#endif              

              ASSERT3 ("Block::refresh_load_extra_face_",
                       "Face if3_er %d %d %d out of bounds",
                       if3_er[0],if3_er[1],if3_er[2],
                       (-1 <= if3_er[0] && if3_er[0] <= 1) &&
                       (-1 <= if3_er[1] && if3_er[1] <= 1) &&
                       (-1 <= if3_er[2] && if3_er[2] <= 1));

            } // level_extra == level
          } // overlap
        } // ! match
      } // while (it_extra.next())
      
    } else if (l_recv) {

      // RECEIVER LOOP OVER EXTRA BLOCKS

      int ef3[3];
      while (it_extra.next(ef3)) {
        
#ifdef DEBUG_ENZO_PROLONG
        CkPrintf ("DEBUG_EXTRA %s recv loop\n",
                  name().c_str());
#endif        
        const Index index_extra = it_extra.index();
        const int   level_extra = it_extra.face_level();
        
        int ec3[3] = {0,0,0};
        if (level_extra > level_face) {
          index_extra.child(level_extra,ec3,ec3+1,ec3+2);
        }

#ifdef DEBUG_ENZO_PROLONG
        CkPrintf ("DEBUG_EXTRA %s recv loop: if3 %d %d %d ic3 %d %d %d\n",
                  name().c_str(),
                  if3[0],if3[1],if3[2],
                  ic3[0],ic3[1],ic3[2]);
        CkPrintf ("DEBUG_EXTRA %s recv loop: ef3 %d %d %d ic3 %d %d %d\n",
                  name().c_str(),
                  ef3[0],ef3[1],ef3[2],
                  ec3[0],ec3[1],ec3[2]);
        CkPrintf ("DEBUG_EXTRA %s recv loop: Ls Lr Le %d %d %d\n",
                  name().c_str(),
                  level_face,level,level_extra);
#endif

        // ... skip extra block if it's the same as the neighbor

        bool l_valid_level = (std::abs(level_extra - level_face) <= 1);

        if (index_extra != index_neighbor && l_valid_level) {

        
#ifdef DEBUG_ENZO_PROLONG
          CkPrintf ("DEBUG_EXTRA %s recv unique\n",name().c_str());
#endif          
          
          // *** count expected receive from Be or be ***

          // ... determine overlap of extra block with intersection region
          const int level_send = level_face;
          int if3_se[3] = {ef3[0]-if3[0],ef3[1]-if3[1],ef3[2]-if3[2] };
          // adjust for fine blocks pointing to neighboring fine block in same parent
          for (int i=0; i<3; i++) {
            if (if3[i] == 0) {
              if (ef3[i] == -1 && ec3[i] == 0) if3_se[i] = 0;
              if (ef3[i] == +1 && ec3[i] == 1) if3_se[i] = 0;
            }
          }
#ifdef DEBUG_ENZO_PROLONG
          CkPrintf ("DEBUG_RECV sr %d %d   %d %d - %d %d\n",
                    if3_se[0],if3_se[1],
                    ef3[0],ef3[1],
                    if3[0],if3[1]);
#endif          
          box_sr.set_block ((level_extra-level_send), if3_se,ec3);
          box_sr.compute_block_start();
        
          int im3[3],ip3[3];
          bool overlap = box_sr.get_limits (im3,ip3,Box::BlockType::extra);
#ifdef DEBUG_ENZO_PROLONG
          box_sr.print("sr overlap recv");
#endif          

#ifdef DEBUG_ENZO_PROLONG
          CkPrintf ("DEBUG_PROLONG recv overlap %d if3 %d %d ef3 %d %d im3 %d %d ip3 %d %d\n",
                    overlap?1:0,
                    if3[0],if3[1],
                    ef3[0],ef3[1],
                    im3[0],im3[1],
                    ip3[0],ip3[1]);
#endif          
          if (overlap) {

            ++count;
            TRACE_COMM_XPCT("extra",name(),refresh.id(),name(index_extra));
            
#ifdef DEBUG_ENZO_PROLONG
            CkPrintf ("DEBUG_EXTRA %s recv overlap count %d\n",
                      name().c_str(),count);
#endif            
        
            if (level_extra == level) {

#ifdef DEBUG_ENZO_PROLONG
              CkPrintf ("DEBUG_EXTRA %s recv extra\n",name().c_str());
#endif              
              // this block receives; extra block is fine

              // handle contribution of this block br to Bs -> be

              // if3 rs
              // ef3 re
              // se = ef3 -
              // int if3_se[3] = { ef3[0]+0.5*ec3[0]-(if3[0]+0.5*ic3[0]),
              //                   ef3[1]+0.5*ec3[1]-(if3[1]+0.5*ic3[1]),
              //                   ef3[2]+0.5*ec3[2]-(if3[2]+0.5*ic3[2])};
            
              // Box br | Bs -> be
              box_se.set_block(+1,if3_se,ec3);
              box_se.set_padding(padding);
              box_se.compute_region();

              int itm3[3],itp3[3];
              box_se.get_limits (itm3,itp3,Box::BlockType::send);

              int if3_sr[3] = { -if3[0], -if3[1], -if3[2] };

              box_se.set_block(0,if3_sr,ic3);
              box_se.compute_block_start();

              int ifm3[3],ifp3[3];
              int iam3[3],iap3[3];
              box_se.get_limits (ifm3,ifp3,Box::BlockType::extra);
              box_se.get_limits (iam3,iap3,Box::BlockType::array);

              int ma3[3] =
                { itp3[0]-itm3[0], itp3[1]-itm3[1], itp3[2]-itm3[2]};

              refresh_extra_send_
                (refresh, index_extra, if3_se,
                 ma3, iam3,iap3, ifm3,ifp3, data()->field());

#ifdef DEBUG_ENZO_PROLONG
              CkPrintf ("DEBUG_PROLONG recv array size  %d %d %d\n",
                        itp3[0]-itm3[0],itp3[1]-itm3[1],itp3[2]-itm3[2]);
              CkPrintf ("DEBUG_PROLONG recv limits ifm3 %d %d %d ifp3 %d %d %d\n",
                        ifm3[0],ifm3[1],ifm3[2],ifp3[0],ifp3[1],ifp3[2]);
              CkPrintf ("DEBUG_PROLONG recv limits iam3 %d %d %d iap3 %d %d %d\n\n",
                        iam3[0],iam3[1],iam3[2],iap3[0],iap3[1],iap3[2]);
#endif              
              
              ASSERT3 ("Block::refresh_load_extra_face_",
                       "Face if3_se %d %d %d out of bounds",
                       if3_se[0],if3_se[1],if3_se[2],
                       ((-1 <= if3_se[0] && if3_se[0] <= 1) &&
                        (-1 <= if3_se[1] && if3_se[1] <= 1) &&
                        (-1 <= if3_se[2] && if3_se[2] <= 1)));

            } // level_extra == level
          } // if (overlap)
        } // if (! match)
      } // while (it_extra.next())
    } // (level > level_face)
  } // (level != level_face)

#ifdef DEBUG_ENZO_PROLONG
  CkPrintf ("DEBUG_PROLONG %s extra count %d\n",name().c_str(),count);
#endif  

  return count;
}

//----------------------------------------------------------------------

void Block::refresh_extra_send_
(Refresh & refresh, Index index_neighbor, int if3[3],
 int m3[3], int iam3[3], int iap3[3], int ifm3[3], int ifp3[3],
 Field field)
{

  MsgRefresh * msg_refresh = new MsgRefresh;
#ifdef TRACE_MSG_REFRESH  
  msg_refresh->set_block_type(name(),"extra");
#endif  

  DataMsg * data_msg = new DataMsg;


  const int id_refresh = refresh.id();
  CHECK_ID(id_refresh);

  data_msg->set_padded_face (if3,m3,iam3,iap3,ifm3, ifp3,
                             refresh.field_list_src(),field);

  msg_refresh->set_refresh_id (id_refresh);
  msg_refresh->set_data_msg (data_msg);

  TRACE_COMM_SEND("extra",this,name(index_neighbor),id_refresh);
  thisProxy[index_neighbor].p_refresh_recv (msg_refresh);
}


//----------------------------------------------------------------------

int Block::refresh_load_particle_faces_ (Refresh & refresh)
{
  const int rank = cello::rank();

  const int npa3[3] = { 4, 4*4, 4*4*4 };
  const int npa = npa3[rank-1];

  ParticleData ** particle_array = new ParticleData*[npa];
  ParticleData ** particle_list = new ParticleData*[npa];
  std::fill_n (particle_array,npa,nullptr);
  std::fill_n (particle_list,npa,nullptr);

  Index * index_list = new Index[npa];
  
  // Sort particles that have left the Block into 4x4x4 array
  // corresponding to neighbors

  int nl = particle_load_faces_
    (npa,particle_list,particle_array, index_list, &refresh);

  // Send particle data to neighbors

  particle_send_(refresh,nl,index_list,particle_list);

  delete [] particle_array;
  delete [] particle_list;
  delete [] index_list;

  return nl;
}

//----------------------------------------------------------------------

void Block::particle_send_
(Refresh & refresh, int nl,Index index_list[], ParticleData * particle_list[])
{

  ParticleDescr * p_descr = cello::particle_descr();

  for (int il=0; il<nl; il++) {

    Index index           = index_list[il];
    ParticleData * p_data = particle_list[il];
    Particle particle_send (p_descr,p_data);
    
    const int id_refresh = refresh.id();
    CHECK_ID(id_refresh);

    ASSERT1 ("Block::particle_send_()",
	     "id_refresh %d of refresh object is out of range",
	     id_refresh,
	     (0 <= id_refresh));

  if (p_data && p_data->num_particles(p_descr)>0) {

      DataMsg * data_msg = new DataMsg;
      data_msg ->set_particle_data(p_data,true);

      MsgRefresh * msg_refresh = new MsgRefresh;
#ifdef TRACE_MSG_REFRESH  
      msg_refresh->set_block_type(name(),"particle");
#endif      
      msg_refresh->set_data_msg (data_msg);
      msg_refresh->set_refresh_id (id_refresh);

      TRACE_COMM_SEND("particle",this,name(index),id_refresh);
      thisProxy[index].p_refresh_recv (msg_refresh);

    } else if (p_data) {
      
      MsgRefresh * msg_refresh = new MsgRefresh;
#ifdef TRACE_MSG_REFRESH  
      msg_refresh->set_block_type(name(),"particle");
#endif      
      msg_refresh->set_data_msg (nullptr);
      msg_refresh->set_refresh_id (id_refresh);

      TRACE_COMM_SEND("particle",this,name(index),id_refresh);
      thisProxy[index].p_refresh_recv (msg_refresh);

      // assert ParticleData object exits but has no particles
      delete p_data;

    }

  }
}

//----------------------------------------------------------------------

int Block::particle_load_faces_ (int npa, 
				 ParticleData * particle_list[],
				 ParticleData * particle_array[],
				 Index index_list[],
				 Refresh *refresh)
{
  // Array elements correspond to child-sized blocks to
  // the left, inside, and right of the main Block.  Particles
  // are assumed to be (well) within this area.
  //
  //     +---+---+---+---+
  //     | 03| 13| 23| 33|
  //     +---+===+===+---+
  //     | 02||  :  || 32|
  //     +---+ - + - +---+
  //     | 01||  :  || 31|
  //     +---+=======+---+
  //     | 00| 10| 20| 30|
  //     +---+---+---+---+
  //
  // Actual neighbors may overlap multiple child-sized blocks.  In
  // that case, we have one ParticleData object per neighbor, but
  // with pointer duplicated.   So if neighbor configuration is:
  //
  //     +---+   5   +---+
  //     | 4 |       | 6 |
  // +---+---+===+===+---+
  // |       ||     ||    
  // |   2   +       +   3
  // |       ||     ||    
  // +-------+=======+-------+
  //         |            
  //     0   |            
  //                 1   
  //
  // Then the particle data array will be:
  //
  //     +---+---+---+---+
  //     | 4 | 5 | 5 | 6 |
  //     +---+===+===+---+
  //     | 2 ||  :  || 3 |
  //     +---+ - + - +---+
  //     | 2 ||  :  || 3 |
  //     +---+=======+---+
  //     | 0 | 1 | 1 | 1 |
  //     +---+---+---+---+

  // ... arrays for updating positions of particles that cross
  // periodic boundaries

  int nl = particle_create_array_neighbors_
    (refresh, particle_array,particle_list,index_list);

  // Scatter particles among particle_data array

  Particle particle (cello::particle_descr(),
		     data()->particle_data());

  std::vector<int> type_list;
  if (refresh->all_particles()) {
    const int nt = particle.num_types();
    type_list.resize(nt);
    for (int i=0; i<nt; i++) type_list[i] = i;
  } else {
    type_list = refresh->particle_list();
  }

  particle_scatter_neighbors_(npa,particle_array,type_list, particle);

  // Update positions particles crossing periodic boundaries

  particle_apply_periodic_update_  (nl,particle_list,refresh);

  return nl;
}

//----------------------------------------------------------------------

int Block::particle_create_array_neighbors_
(Refresh * refresh, 
 ParticleData * particle_array[],
 ParticleData * particle_list[],
 Index index_list[])
{ 
  const int rank = cello::rank();
  const int level = this->level();

  const int min_face_rank = refresh->min_face_rank();

  ItNeighbor it_neighbor =
    this->it_neighbor(min_face_rank,index_, neighbor_leaf,0,0);

  int il = 0;

  int if3[3];
  for (il=0; it_neighbor.next(if3); il++) {

    const int level_face = it_neighbor.face_level();
    TRACE_COMM_XPCT("particle",name(),refresh->id(),name(it_neighbor.index()));

    int ic3[3] = {0,0,0};

    const int refresh_type = 
      (level_face == level - 1) ? refresh_coarse :
      (level_face == level)     ? refresh_same :
      (level_face == level + 1) ? refresh_fine : refresh_unknown;

    if (refresh_type==refresh_coarse) {
      // coarse neighbor: need index of self in parent
      index_.child(index_.level(),ic3,ic3+1,ic3+2);
    } else if (refresh_type==refresh_fine) {
      // fine neighbor: need index of child in self
      it_neighbor.child(ic3);
    }
    // (else same-level neighbor: don't need child)

    int index_lower[3] = {0,0,0};
    int index_upper[3] = {1,1,1};
    refresh->get_particle_bin_limits
      (rank,refresh_type,if3,ic3,index_lower,index_upper);

    ParticleData * pd = new ParticleData;

    ParticleDescr * p_descr = cello::particle_descr();

    pd->allocate(p_descr);

    particle_list[il] = pd;

    index_list[il] = it_neighbor.index();

    for (int iz=index_lower[2]; iz<index_upper[2]; iz++) {
      for (int iy=index_lower[1]; iy<index_upper[1]; iy++) {
	for (int ix=index_lower[0]; ix<index_upper[0]; ix++) {
	  int i=ix + 4*(iy + 4*iz);
	  particle_array[i] = pd;
	}
      }
    }
  }
  
  return il;
}

//----------------------------------------------------------------------

void Block::particle_determine_periodic_update_
(int * index_lower, int * index_upper,
 double *dpx, double *dpy, double *dpz)
{
  //     ... domain extents
  double dxm,dym,dzm;
  double dxp,dyp,dzp;

  cello::hierarchy()->lower(&dxm,&dym,&dzm);
  cello::hierarchy()->upper(&dxp,&dyp,&dzp);

  //     ... periodicity
  bool p3[3];
  periodicity(p3);

  //     ... boundary
  bool b32[3][2];
  is_on_boundary (b32);

  const int rank = cello::rank();

  // Update (dpx,dpy,dpz) position correction if periodic domain
  // boundary is crossed

  if (rank >= 1) {
    if (index_lower[0]==0 && b32[0][0] && p3[0]) (*dpx) = +(dxp - dxm);
    if (index_upper[0]==4 && b32[0][1] && p3[0]) (*dpx) = -(dxp - dxm);
  }
  if (rank >= 2) {
    if (index_lower[1]==0 && b32[1][0] && p3[1]) (*dpy) = +(dyp - dym);
    if (index_upper[1]==4 && b32[1][1] && p3[1]) (*dpy) = -(dyp - dym);
  }
  if (rank >= 3) {
    if (index_lower[2]==0 && b32[2][0] && p3[2]) (*dpz) = +(dzp - dzm);
    if (index_upper[2]==4 && b32[2][1] && p3[2]) (*dpz) = -(dzp - dzm);
  }
}

//----------------------------------------------------------------------

void Block::particle_apply_periodic_update_
(int nl, ParticleData * particle_list[], Refresh * refresh)
{

  const int rank = cello::rank();
  const int level = this->level();
  const int min_face_rank = refresh->min_face_rank();

  std::vector<double> dpx(nl,0.0);
  std::vector<double> dpy(nl,0.0);
  std::vector<double> dpz(nl,0.0);

  // Compute position updates for particles crossing periodic boundaries

  ItNeighbor it_neighbor =
    this->it_neighbor(min_face_rank,index_, neighbor_leaf,0,0);

  int il=0;

  int if3[3];
  while (it_neighbor.next(if3)) {

    const int level_face = it_neighbor.face_level();

    int ic3[3];
    it_neighbor.child(ic3);

    const int refresh_type = 
      (level_face == level - 1) ? refresh_coarse :
      (level_face == level)     ? refresh_same :
      (level_face == level + 1) ? refresh_fine : refresh_unknown;

    int index_lower[3] = {0,0,0};
    int index_upper[3] = {1,1,1};
    refresh->get_particle_bin_limits
      (rank,refresh_type,if3,ic3,index_lower,index_upper);

    // ASSERT: il < nl
    particle_determine_periodic_update_
      (index_lower,index_upper,&dpx[il],&dpy[il],&dpz[il]);

    il++;

  }

  ParticleDescr * p_descr = cello::particle_descr();

  // Apply the updates to the list of particles

  for (int il=0; il<nl; il++) {

    ParticleData * p_data = particle_list[il];
    Particle particle_neighbor (p_descr,p_data);

    if ( ((rank >= 1) && dpx[il] != 0.0) ||
	 ((rank >= 2) && dpy[il] != 0.0) ||
	 ((rank >= 3) && dpz[il] != 0.0) ) {
	
      // ... for each particle type
      const int nt = particle_neighbor.num_types();
      for (int it=0; it<nt; it++) {

	// ... for each batch of particles
	const int nb = particle_neighbor.num_batches(it);
	for (int ib=0; ib<nb; ib++) {

	  particle_neighbor.position_update (it,ib,dpx[il],dpy[il],dpz[il]);

	}
      }
    }
  }
}
//----------------------------------------------------------------------

void Block::particle_scatter_neighbors_
(int npa,
 ParticleData * particle_array[],
 std::vector<int> & type_list,
 Particle particle)
{
  const int rank = cello::rank();

  //     ... get Block bounds 
  double xm,ym,zm;
  double xp,yp,zp;
  lower(&xm,&ym,&zm);
  upper(&xp,&yp,&zp);

  // find block center (x0,y0,z0) and width (xl,yl,zl)
  const double x0 = 0.5*(xm+xp);
  const double y0 = 0.5*(ym+yp);
  const double z0 = 0.5*(zm+zp);
  const double xl = xp-xm;
  const double yl = yp-ym;
  const double zl = zp-zm;

  int count = 0;
  // ...for each particle type to be moved

  for (auto it_type=type_list.begin(); it_type!=type_list.end(); it_type++) {

    int it = *it_type;

    const int ia_x  = particle.attribute_position(it,0);

    // (...positions may use absolute coordinates (float) or
    // block-local coordinates (int))
    const bool is_float = 
      (cello::type_is_float(particle.attribute_type(it,ia_x)));

    // (...stride may be != 1 if particle attributes are interleaved)
    const int d  = particle.stride(it,ia_x);

    // ...for each batch of particles

    const int nb = particle.num_batches(it);

    for (int ib=0; ib<nb; ib++) {

      const int np = particle.num_particles(it,ib);

      // ...extract particle position arrays

      std::vector<double> xa(np,0.0);
      std::vector<double> ya(np,0.0);
      std::vector<double> za(np,0.0);

      particle.position(it,ib,xa.data(),ya.data(),za.data());

      // ...initialize mask used for scatter and delete
      // ...and corresponding particle indices

      bool * mask = new bool[np];
      int * index = new int[np];
      
      for (int ip=0; ip<np; ip++) {

	double x = is_float ? 2.0*(xa[ip*d]-x0)/xl : xa[ip*d];
	double y = is_float ? 2.0*(ya[ip*d]-y0)/yl : ya[ip*d];
	double z = is_float ? 2.0*(za[ip*d]-z0)/zl : za[ip*d];

	int ix = (rank >= 1) ? (x + 2) : 0;
	int iy = (rank >= 2) ? (y + 2) : 0;
	int iz = (rank >= 3) ? (z + 2) : 0;

	if (! (0 <= ix && ix < 4) ||
	    ! (0 <= iy && iy < 4) ||
	    ! (0 <= iz && iz < 4)) {
	  
	  CkPrintf ("%d ix iy iz %d %d %d\n",CkMyPe(),ix,iy,iz);
	  CkPrintf ("%d x y z %f %f %f\n",CkMyPe(),x,y,z);
	  CkPrintf ("%d xa ya za %f %f %f\n",CkMyPe(),xa[ip*d],ya[ip*d],za[ip*d]);
	  CkPrintf ("%d xm ym zm %f %f %f\n",CkMyPe(),xm,ym,zm);
	  CkPrintf ("%d xp yp zp %f %f %f\n",CkMyPe(),xp,yp,zp);
	  ERROR3 ("Block::particle_scatter_neighbors_",
		  "particle indices (ix,iy,iz) = (%d,%d,%d) out of bounds",
		  ix,iy,iz);
	}

	const int i = ix + 4*(iy + 4*iz);
	index[ip] = i;
	bool in_block = true;
	in_block = in_block && (!(rank >= 1) || (1 <= ix && ix <= 2));
	in_block = in_block && (!(rank >= 2) || (1 <= iy && iy <= 2));
	in_block = in_block && (!(rank >= 3) || (1 <= iz && iz <= 2));
	mask[ip] = ! in_block;
      }

      // ...scatter particles to particle array
      particle.scatter (it,ib, np, mask, index, npa, particle_array);
      // ... delete scattered particles
      count += particle.delete_particles (it,ib,mask);

      delete [] mask;
      delete [] index;
    }
  }

  cello::simulation()->data_delete_particles(count);

}

//----------------------------------------------------------------------

int Block::refresh_load_flux_faces_ (Refresh & refresh)
{
  int count = 0;

  const int min_face_rank = cello::rank() - 1;
  const int neighbor_type = neighbor_leaf;
  
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
      (level_face < level) ? refresh_coarse :
      (level_face > level) ? refresh_fine : refresh_same;

    refresh_load_flux_face_
      (refresh,refresh_type,index_neighbor,if3,ic3);

    TRACE_COMM_XPCT("fluxes",name(),refresh.id(),name(index_neighbor));
    ++count;

  }

  return count;
}

//----------------------------------------------------------------------

void Block::refresh_load_flux_face_
( Refresh & refresh,
  int refresh_type,
  Index index_neighbor,
  int if3[3],
  int ic3[3])

{
  // ... coarse neighbor requires child index of self in parent
  if (refresh_type == refresh_coarse) {
    index_.child(index_.level(),ic3,ic3+1,ic3+2);
  }

  // ... copy field ghosts to array using FieldFace object

  const int axis = (if3[0]!=0) ? 0 : (if3[1]!=0) ? 1 : 2;
  const int face = (if3[axis]==-1) ? 0 : 1;


  // neighbor is coarser
  DataMsg * data_msg = new DataMsg;
  FluxData * flux_data = data()->flux_data();
  
  const bool is_new = true;
  if (refresh_type == refresh_coarse) {
    // neighbor is coarser
    const int nf = flux_data->num_fields();
    data_msg -> set_num_face_fluxes(nf);
    for (int i=0; i<nf; i++) {
      FaceFluxes * face_fluxes = new FaceFluxes
        (*flux_data->block_fluxes(axis,face,i));
      face_fluxes->coarsen(ic3[0],ic3[1],ic3[2],cello::rank());
      data_msg -> set_face_fluxes (i,face_fluxes, is_new);
    }
  } else {
    data_msg -> set_num_face_fluxes(0);
  }

  const int id_refresh = refresh.id();
  CHECK_ID(id_refresh);
  
  ASSERT1 ("Block::refresh_load_flux_face_()",
           "id_refresh %d of refresh object is out of range",
           id_refresh,
           (0 <= id_refresh));

  MsgRefresh * msg_refresh = new MsgRefresh;
#ifdef TRACE_MSG_REFRESH  
  msg_refresh->set_block_type(name(),"fluxes");
#endif  
  msg_refresh->set_data_msg (data_msg);
  msg_refresh->set_refresh_id (id_refresh);

  TRACE_COMM_SEND("fluxes",this,name(index_neighbor),id_refresh);
  thisProxy[index_neighbor].p_refresh_recv (msg_refresh);

}
