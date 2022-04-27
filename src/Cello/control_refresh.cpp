// See LICENSE_CELLO file for license and copyright information

/// @file     control_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2019-05-23
/// @brief    Charm-related functions associated with refreshing ghost zones
/// @ingroup  Control

// #define TRACE_LOAD_FACE
// #define TRACE_PROLONG

// #define PRINT_COARSE_FIELD
// #define DEBUG_ARRAY
// #define DEBUG_ARRAY_CYCLE 0
// #define DEBUG_PRINT true
// #define DEBUG_PRINT_BLOCK "B1:0_1:1"
// #define DEBUG_BOX

#include "simulation.hpp"
#include "mesh.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

#define CHECK_ID(ID) ASSERT1 ("CHECK_ID","Invalid id %d",ID,(ID>=0));

#ifdef TRACE_PROLONG
#  undef TRACE_PROLONG
#  define TRACE_PROLONG(MSG,PROLONG,mf3,if3,nf3,mc3,ic3,nc3)            \
  CkPrintf ("TRACE_PROLONG %s:%d %s %s mf %d %d %d nf %d %d %d if %d %d %d\n",__FILE__,__LINE__,MSG, PROLONG->name().c_str(), \
            mf3[0],mf3[1],mf3[2],nf3[0],nf3[1],nf3[2],if3[0],if3[1],if3[2]); \
  CkPrintf ("TRACE_PROLONG %s:%d %s %s mc %d %d %d nc %d %d %d ic %d %d %d\n",__FILE__,__LINE__,MSG, PROLONG->name().c_str(), \
    mc3[0],mc3[1],mc3[2],nc3[0],nc3[1],nc3[2],ic3[0],ic3[1],ic3[2]);    \
  
#else
#  undef TRACE_PROLONG
#  define TRACE_PROLONG(MSG,PROLONG,mf3,if3,nf3,mc3,ic3,nc3) /* ... */
#endif

//======================================================================

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
      count_particle = refresh_load_particle_faces_(*refresh,
						    refresh->particles_are_copied());
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

  for (size_t id_msg=0;
       id_msg<refresh_msg_list_[id_refresh].size();
       id_msg++) {

    MsgRefresh * msg = refresh_msg_list_[id_refresh][id_msg];

    // unpack message data into Block data
    msg->update(data());

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

    // complete any outstanding data operations (e.g. multi-block
    // interpolations)

    refresh_coarse_apply_(refresh);

    // Call callback

    refresh_exit(*refresh);
  }
}

//----------------------------------------------------------------------

void Block::p_refresh_recv (MsgRefresh * msg_refresh)
{
  const int id_refresh = msg_refresh->id_refresh();
  CHECK_ID(id_refresh);

  Sync * sync = sync_(id_refresh);

  if (sync->state() == RefreshState::READY) {

    // unpack message data into Block data if ready
    msg_refresh->update(data());

    delete msg_refresh;

    sync->advance();

    // check if it's the last message processed
    refresh_check_done(id_refresh);

  } else {

    // save message if not ready
    refresh_msg_list_[id_refresh].push_back(msg_refresh);

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
      this->it_neighbor(index_,min_face_rank,
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
      Prolong * prolong = refresh.prolong();
      int pad = refresh.coarse_padding(prolong);

      if (pad == 0) {
        refresh_load_field_face_
          (refresh,refresh_type,index_neighbor,if3,ic3);
        ++count;
      } else {
        if (level_face == level) {
          refresh_load_field_face_
            (refresh,refresh_type,index_neighbor,if3,ic3);
          ++count;
        } else {
          count += refresh_load_coarse_face_
            (refresh,refresh_type,index_neighbor,if3,ic3);
        }
        if (level_face < level) {
          refresh_load_field_face_
            (refresh,refresh_type,index_neighbor,if3,ic3);
        } else if (level_face > level) {
          count ++;
        }

      }
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

      }

    }
  }
  return count;
}

//----------------------------------------------------------------------

void Block::refresh_load_field_face_
( Refresh & refresh,  int refresh_type,
  Index index_neighbor,  int if3[3], int ic3[3])
{
  // create refresh message

  MsgRefresh * msg_refresh = new MsgRefresh;

  // create field face
  if (refresh_type == refresh_coarse) {
    index_.child(index_.level(),ic3,ic3+1,ic3+2);
  }
  int g3[3] = {0,0,0};
  FieldFace * field_face = create_face
    (if3, ic3, g3, refresh_type, &refresh,false);

  // create data message
  DataMsg * data_msg = new DataMsg;
  // initialize data message
  data_msg -> set_field_face (field_face,true);
  data_msg -> set_field_data (data()->field_data(),false);

  // initialize refresh message
  msg_refresh->set_refresh_id (refresh.id());
  msg_refresh->set_data_msg (data_msg);

  thisProxy[index_neighbor].p_refresh_recv (msg_refresh);

}

//----------------------------------------------------------------------

int Block::refresh_load_coarse_face_
(Refresh refresh, int refresh_type,
 Index index_neighbor, int if3[3], int ic3[3])
{
  const int level_face = index_neighbor.level();
  const int level = index_.level();
  int count = 0;

  Prolong * prolong = refresh.prolong();
  const int pad = refresh.coarse_padding(prolong);

  if ((pad > 0) && (level != level_face)) {

    // Create box_face

    const int rank = cello::rank();
    int n3[3];
    data()->field().size(n3,n3+1,n3+2);
    const int g = refresh.ghost_depth();
    int g3[3] = {(rank >= 1) ? g : 0,
                 (rank >= 2) ? g : 0,
                 (rank >= 3) ? g : 0};

    // Boxes used Bs -> br, Be -> br, Bs -> be
    // (s send, r receive, e extra)
    Box box_sr (rank,n3,g3);
    Box box_se (rank,n3,g3);
    Box box_er (rank,n3,g3);

    // Create iterator over extra blocks

    ItNeighbor it_extra =
      this->it_neighbor(index_,refresh.min_face_rank(),
                        refresh.neighbor_type(),
                        cello::config()->mesh_min_level,
                        refresh.root_level());

      // ... determine intersection region

    const bool l_send = (level < level_face);
    const bool l_recv = (level > level_face);

    int jf3[3] = { l_send ? if3[0] : -if3[0],
                   l_send ? if3[1] : -if3[1],
                   l_send ? if3[2] : -if3[2] };

    box_sr.set_block(BoxType_receive,+1,jf3,ic3);
    box_sr.set_padding(pad);

    // ... adjust send-ghost depth for accumulate

    refresh.box_accumulate_adjust(&box_sr,if3,g3);

    box_sr.compute_region();

    // SEND MAIN SEND->RECV VIA COARSE ARRAY

    int iam3[3],iap3[3];
    int ifms3[3],ifps3[3];
    int ifmr3[3],ifpr3[3];

    if (l_send) {

      bool lpad;

      box_sr.get_start_stop
        (iam3,iap3,BlockType::extra,BlockType::receive_coarse,lpad=true);
      box_sr.get_start_stop
        (ifms3,ifps3,BlockType::extra,BlockType::send,lpad=true);
      box_sr.get_start_stop
        (ifmr3,ifpr3,BlockType::extra,BlockType::receive,lpad=false);

      refresh_coarse_send_
        (index_neighbor,
         data()->field(), refresh,
         iam3,iap3,ifms3,ifps3,ifmr3,ifpr3,"SR");

    } else if (l_recv) {

      // only count receives
      count ++;

    }

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

        // ... skip extra block if it's the same as the neighbor

        bool l_valid_level = (std::abs(level_extra - level) <= 1);

        if (index_extra != index_neighbor && l_valid_level) {

          // ... determine overlap of extra block with intersection region
          const int level_send = level;
          box_sr.set_block (BoxType_extra,(level_extra-level_send), ef3,ec3);
          int tm3[3],tp3[3];
          bool lpad;

          bool overlap = box_sr.get_start_stop
            (tm3,tp3,BlockType::extra,BlockType::extra,lpad=true);
          if (overlap) {

          if (level_extra == level) {

              // this block sends; extra block is coarse

              // handle contribution of this block Bs to Be -> br

              int if3_er[3] = { if3[0]-ef3[0], if3[1]-ef3[1], if3[2]-ef3[2] };


              // Box Bs | Be -> br
              box_er.set_block(BoxType_receive,+1,if3_er,ic3);
              box_er.set_padding(pad);

              // ... adjust send-ghost depth for accumulate
              refresh.box_accumulate_adjust(&box_er,if3_er,g3);

              box_er.compute_region();


              int if3_es[3] = {-ef3[0], -ef3[1], -ef3[2] };

              box_er.set_block(BoxType_extra,0,if3_es,ic3); // ic3 ignored

              bool lpad;

              box_er.get_start_stop
                (iam3,iap3,BlockType::extra,BlockType::receive_coarse,lpad=true);
              box_er.get_start_stop
                (ifms3,ifps3,BlockType::extra,BlockType::extra,lpad=true);
              box_er.get_start_stop
                (ifmr3,ifpr3,BlockType::extra,BlockType::receive,lpad=false);

              refresh_coarse_send_
                (index_neighbor,
                 data()->field(), refresh,
                 iam3,iap3,ifms3,ifps3,ifmr3,ifpr3,"ER");

              ASSERT3 ("Block::refresh_load_coarse_face_",
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

        const Index index_extra = it_extra.index();
        const int   level_extra = it_extra.face_level();

        int ec3[3] = {0,0,0};
        if (level_extra > level_face) {
          index_extra.child(level_extra,ec3,ec3+1,ec3+2);
        }

        // ... skip extra block if it's the same as the neighbor

        bool l_valid_level = (std::abs(level_extra - level_face) <= 1);

        if (index_extra != index_neighbor && l_valid_level) {


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
          box_sr.set_block (BoxType_extra,(level_extra-level_send), if3_se,ec3);

          int tm3[3],tp3[3];
          bool lpad;

          bool overlap = box_sr.get_start_stop
            (tm3,tp3,BlockType::extra,BlockType::extra,lpad = true);
          if (overlap) {

            ++count;

            if (level_extra == level) {

              // handle contribution of this block br to Bs -> be
              // Box br | Bs -> be
              box_se.set_block(BoxType_receive,+1,if3_se,ec3);
              box_se.set_padding(pad);

              // ... adjust send-ghost depth for accumulate
              refresh.box_accumulate_adjust(&box_se,if3_se,g3);

              box_se.compute_region();

              int if3_sr[3] = { -if3[0], -if3[1], -if3[2] };
              box_se.set_block(BoxType_extra,+1,if3_sr,ic3);

              int iam3[3],iap3[3];
              int ifms3[3],ifps3[3];
              int ifmr3[3],ifpr3[3];
              bool lpad;

              box_se.get_start_stop
                (iam3,iap3,BlockType::extra,BlockType::receive_coarse,lpad=true);
              box_se.get_start_stop
                (ifms3,ifps3,BlockType::extra,BlockType::extra,lpad=true);
              box_se.get_start_stop
                (ifmr3,ifpr3,BlockType::extra,BlockType::receive,lpad=false);

              int ma3[3];
              box_se.get_region_size(ma3);

              refresh_coarse_send_
                (index_extra,
                 data()->field(),refresh,
                 iam3,iap3, ifms3,ifps3, ifmr3,ifpr3,"SE");

              ASSERT3 ("Block::refresh_load_coarse_face_",
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

  return count;
}

//----------------------------------------------------------------------


int Block::delete_non_local_particles_(int it){

  Particle particle (cello::particle_descr(),
		     data()->particle_data());

  ASSERT1("Block::clean_up_particles_",
	  "This function has been called for particle type"
	  " %s, which has no is_copy attribute",
	  particle.type_name(it),
	  particle.is_attribute(it,"is_copy"));

  const int rank = cello::rank();

  // Get Block bounds
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

  const int ia_x  = particle.attribute_position(it,0);

  // (...positions may use absolute coordinates (float) or
  // block-local coordinates (int))
  const bool is_float =
    (cello::type_is_float(particle.attribute_type(it,ia_x)));
  
  // (...stride may be != 1 if particle attributes are interleaved)
  const int d  = particle.stride(it,ia_x);

  // Need pointer to is_copy attribute
  const int ia_copy = particle.attribute_index(it,"is_copy");
  const int d_copy   = particle.stride(it, ia_copy);
  int64_t * is_copy=0;

  // Loop over batches
  const int nb = particle.num_batches(it);
  for (int ib = 0; ib<nb; ib++){
    const int np = particle.num_particles(it,ib);

    if (np == 0) continue;

    // ...extract particle position arrays
	
    std::vector<double> xa(np,0.0);
    std::vector<double> ya(np,0.0);
    std::vector<double> za(np,0.0);

    particle.position(it,ib,xa.data(),ya.data(),za.data());
	
    is_copy = (int64_t *) particle.attribute_array(it, ia_copy, ib);

    // Mask will be used to delete particles
    bool * mask = new bool[np];
    for( int ip=0; ip<np; ip++){

      // We check if particles are within the bounds of the block
      // look at block scatter children for help?
      double x = is_float ? 2.0*(xa[ip*d]-x0)/xl : xa[ip*d];
      double y = is_float ? 2.0*(ya[ip*d]-y0)/yl : ya[ip*d];
      double z = is_float ? 2.0*(za[ip*d]-z0)/zl : za[ip*d];
      
      int ix = (rank >= 1) ? (x + 2) : 0;
      int iy = (rank >= 2) ? (y + 2) : 0;
      int iz = (rank >= 3) ? (z + 2) : 0;

      bool in_block = true;

      // If particle is not in block, it is tagged for deletion
      in_block = in_block && (!(rank >= 1) || (1 <= ix && ix <= 2));
      in_block = in_block && (!(rank >= 2) || (1 <= iy && iy <= 2));
      in_block = in_block && (!(rank >= 3) || (1 <= iz && iz <= 2));
      mask[ip] = ! in_block;

      // Also, we set is_copy = false for particles left behind, for which
      // in_block = true.
      is_copy[ip*d_copy] = !in_block;
    }

    count += particle.delete_particles(it,ib,mask);

    delete [] mask;
  }

  return count;
}

//----------------------------------------------------------------------

void Block::refresh_coarse_send_
(Index index_neighbor,
 Field field,Refresh & refresh,
 int iam3[3], int iap3[3],
 int ifms3[3], int ifps3[3],
 int ifmr3[3], int ifpr3[3],
 std::string debug)
{
  MsgRefresh * msg_refresh = new MsgRefresh;

  DataMsg * data_msg = new DataMsg;

  const int id_refresh = refresh.id();
  CHECK_ID(id_refresh);
  data_msg->set_coarse_array
    (field, iam3,iap3,ifms3,ifps3,ifmr3,ifpr3,
     refresh.field_list_src(),
     refresh.field_list_dst(),
     name(index_neighbor));

  msg_refresh->set_refresh_id (id_refresh);
  msg_refresh->set_data_msg (data_msg);

  thisProxy[index_neighbor].p_refresh_recv (msg_refresh);
}

//----------------------------------------------------------------------

void Block::refresh_coarse_apply_ (Refresh * refresh)
{
  Prolong * prolong = refresh->prolong();

  const int pad = refresh->coarse_padding(prolong);

  if (pad > 0) {

    const int min_face_rank = refresh->min_face_rank();
    const int neighbor_type = refresh->neighbor_type();
    const int root_level    = refresh->root_level();
    const int min_level     = cello::config()->mesh_min_level;

    if (neighbor_type == neighbor_leaf ||
        neighbor_type == neighbor_tree) {

      ItNeighbor it_neighbor =
        this->it_neighbor(index_,min_face_rank,neighbor_type,
                          min_level,root_level);

      const int level = this->level();

      int if3[3];
      while (it_neighbor.next(if3)) {

        int of3[3] = {-if3[0],-if3[1],-if3[2]};

        const int level_face = it_neighbor.face_level();

        if (level == level_face + 1) {

          const auto & field_list_src = refresh->field_list_src();
          const auto & field_list_dst = refresh->field_list_dst();
          const int nf = field_list_src.size();

          // Create Box to find loop limits for copying own values
          int n3[3];
          Field field = this->data()->field();
          field.size(n3,n3+1,n3+2);
          const int g = refresh->ghost_depth();
          const int rank = cello::rank();
          int g3[3] = {(rank >= 1) ? g : 0,
                       (rank >= 2) ? g : 0,
                       (rank >= 3) ? g : 0};
          Box box_r (rank,n3,g3);
          int ic3[3];
          index_.child(level,ic3,ic3+1,ic3+2);
          box_r.set_block (BoxType_receive,+1, of3,ic3);
          box_r.set_padding(pad);
          int gr3[3] = {0,0,0};
          box_r.set_recv_ghosts(gr3); // non-ghost block values only
          refresh->box_accumulate_adjust(&box_r,of3,g3);
          box_r.compute_region();

          for (int i_f=0; i_f<nf; i_f++) {

            // ... adjust send-ghost depth for accumulate
            int i3_c[3], n3_c[3];
            int i3_f[3], n3_f[3];
            bool lpad;

            box_r.get_start_size
              (i3_f,n3_f,BlockType::receive,BlockType::receive,lpad=false);
            box_r.get_start_size
              (i3_c,n3_c,BlockType::receive,BlockType::receive_coarse,lpad=true);
            const int index_field_src = field_list_src[i_f];
            const int index_field_dst = field_list_dst[i_f];

            int m3_c[3];
            int m3_f[3];
            field.coarse_dimensions(i_f,m3_c,m3_c+1,m3_c+2);
            field.dimensions(index_field_src,&m3_f[0],&m3_f[1],&m3_f[2]);

            cello_float * field_values_src =
              (cello_float *) field.values(index_field_src);
            cello_float * field_values_dst =
              (cello_float *) field.values(index_field_dst);
            cello_float * coarse_field_src =
              (cello_float *) field.coarse_values(index_field_src);

            const float r = n3_f[0] / n3_c[0];
            const cello_float rr = (r == 1) ? 1.0 : 1.0/cello::num_children();

#ifdef CHECK
            ASSERT1 ("Block::refresh_coarse_apply",
                     "Field-to-coarse array axis ratio r=%g is not 1.0 or 2.0",
                     r, (r==1.0 || r==2.0));
#endif

            const int if0 = i3_f[0] + m3_f[0]*(i3_f[1] + m3_f[1]*i3_f[2]);
            const int ic0 = i3_c[0] + m3_c[0]*(i3_c[1] + m3_c[1]*i3_c[2]);

            // clear first to avoid double-counting overlapped subregions
            for (int kz=0; kz<n3_f[2]; kz++) {
              const int kzr=kz/r;
              for (int ky=0; ky<n3_f[1]; ky++) {
                const int kyr=ky/r;
                for (int kx=0; kx<n3_f[0]; kx++) {
                  const int kxr=kx/r;
                  int kc = ic0 + kxr + m3_c[0]*(kyr + m3_c[1]*kzr);
                  coarse_field_src[kc] = 0;
                }
              }
            }
            for (int kz=0; kz<n3_f[2]; kz++) {
              const int kzr=kz/r;
              for (int ky=0; ky<n3_f[1]; ky++) {
                const int kyr=ky/r;
                for (int kx=0; kx<n3_f[0]; kx++) {
                  const int kxr=kx/r;
                  int kc = ic0 + kxr + m3_c[0]*(kyr + m3_c[1]*kzr);
                  int kf = if0 + kx + m3_f[0]*(ky + m3_f[1]*kz);
                  coarse_field_src[kc] += rr*field_values_src[kf];
                }
              }
            }

            int ip3_c[3],np3_c[3];
            int ip3_f[3],np3_f[3];

            Box box_p (rank,n3,g3);
            int ic3[3];
            index_.child(level,ic3,ic3+1,ic3+2);
            box_p.set_block (BoxType_receive,+1, of3,ic3);
            box_p.set_padding(pad);
            box_p.set_recv_ghosts(g3); // reset recv ghosts to default
            refresh->box_accumulate_adjust(&box_p,of3,g3);

            box_p.compute_region();

            box_p.get_start_size
              (ip3_c,np3_c,BlockType::none,BlockType::receive_coarse,lpad=true);
            box_p.get_start_size
              (ip3_f,np3_f,BlockType::none,BlockType::receive,lpad=false);

            prolong->array_sizes_valid (n3_f,n3_c);
            const bool accumulate = refresh->accumulate(i_f);
            TRACE_PROLONG("coarse_apply",prolong, m3_f,ip3_f,np3_f, m3_c,ip3_c,np3_c);
            prolong->apply(default_precision,
                           field_values_dst, m3_f, ip3_f, np3_f,
                           coarse_field_src, m3_c, ip3_c, np3_c,
                           accumulate);

            if (accumulate) {
              DEBUG_PRINT_ARRAY0("refresh_coarse_apply coarse_field",coarse_field,m3_c,np3_c,ip3_c);
              DEBUG_PRINT_ARRAY0("refresh_field_apply field_values",field_values_dst,m3_f,np3_f,ip3_f);
            }
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------

int Block::refresh_load_particle_faces_ (Refresh & refresh, const bool copy)
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
    (npa,particle_list,particle_array, index_list, &refresh, copy);

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
      msg_refresh->set_data_msg (data_msg);
      msg_refresh->set_refresh_id (id_refresh);

      thisProxy[index].p_refresh_recv (msg_refresh);

    } else if (p_data) {

      MsgRefresh * msg_refresh = new MsgRefresh;
      msg_refresh->set_data_msg (nullptr);
      msg_refresh->set_refresh_id (id_refresh);

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
				 Refresh *refresh,
                                 const bool copy)
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

  particle_scatter_neighbors_(npa,particle_array,type_list, particle, copy);

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
    this->it_neighbor(index_, min_face_rank,neighbor_leaf,0,0);

  int il = 0;

  int if3[3];
  for (il=0; it_neighbor.next(if3); il++) {

    const int level_face = it_neighbor.face_level();

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
  int p3[3];
  //  cello::hierarchy()->get_periodicity(p3);
//  int p3[3];
  cello::hierarchy()->get_periodicity(p3,p3+1,p3+2);

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
    this->it_neighbor(index_, min_face_rank,neighbor_leaf,0,0);

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
 Particle particle,
 const bool copy)
{
  
  if (copy){

     // Loop over particle types
     for (auto it_type=type_list.begin(); it_type!=type_list.end(); it_type++) {

       int it = *it_type;
       
       ASSERT1("Block::particle_scatter_neighbors_",
	       "Trying to copy particle type %s, but it has no"
	       "is_copy attribute",
	       particle.type_name(it),
	       particle.is_attribute(it,"is_copy"));

       int ia_copy = particle.attribute_index(it, "is_copy");
       int d_copy = particle.stride(it,ia_copy);
       int64_t * is_copy;
       
       const int nb = particle.num_batches(it);

       // Loop over batches
       for (int ib=0; ib<nb; ib++) {

	 const int np = particle.num_particles(it,ib);

	 if (np == 0) continue;

	 is_copy = (int64_t *) particle.attribute_array(it, ia_copy, ib);

	 // ...initialize mask used for copying
	 bool * mask = new bool[np];
	 // Index array not needed for copying
	 int * index = nullptr;

	 // Loop over particles in this batch and fill in the mask
	 for (int ip=0; ip<np; ip++) mask[ip] = !is_copy[ip*d_copy];
	 
	 // ...scatter particles to particle array
	 particle.scatter  (it,ib,np,mask,index,npa,particle_array, copy);
	 
	 delete [] mask;
       } // Loop over batches
     } // Loop over particle types
  } // if (copy)

  else {
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
	
	if (np == 0) continue;
	
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

	  // look at block scatter children for help?
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
	particle.scatter  (it,ib,np,mask,index,npa,particle_array, copy);
	
	// ... delete scattered particles if moved
	count += particle.delete_particles (it,ib,mask);
	
	delete [] mask;
	delete [] index;
      } // Loop over batches
    } // Loop over particle types
    
    cello::simulation()->data_delete_particles(count);

  } // else
  return;
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
    this->it_neighbor(index_,min_face_rank,
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
  msg_refresh->set_data_msg (data_msg);
  msg_refresh->set_refresh_id (id_refresh);

  thisProxy[index_neighbor].p_refresh_recv (msg_refresh);

}
