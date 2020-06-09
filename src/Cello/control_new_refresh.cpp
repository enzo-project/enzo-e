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

#define CHECK_ID(ID)				\
  ASSERT1 ("CHECK_ID",				\
	   "Invalid id %d",ID,			\
	   ID>=0);				\

void Block::new_refresh_start (int id_refresh, int callback)
{
  CHECK_ID(id_refresh);

  RefreshState & state = new_refresh_state_list_[id_refresh];
  Refresh * refresh    = cello::refresh(id_refresh);

  // Send field and/or particle data associated with the given refresh
  // object to corresponding neighbors

  if ( refresh->is_active() ) {

    ASSERT1 ("Block::new_refresh_start()",
	     "refresh[%d] state is not inactive",
	     id_refresh,
	     (state == RefreshState::INACTIVE));

    state = RefreshState::ACTIVE;

    int count = 0;

    // send Field face data

    if (refresh->any_fields()) {
      count += new_refresh_load_field_faces_ (*refresh);
    }

    // send Particle face data
    if (refresh->any_particles()){
      count += new_refresh_load_particle_faces_(*refresh);
    }

    if (refresh->any_particles_copy()){
      new_refresh_delete_particle_copies_(refresh);
      count += new_refresh_load_particle_copy_(*refresh);
    }

    Sync & sync = new_refresh_sync_list_[id_refresh];

    // Make sure sync counter is not active
    ASSERT4 ("Block::new_refresh_start()",
	     "refresh[%d] sync object %p is active (%d/%d)",
	     id_refresh, &sync, sync.value(), sync.stop(),
	     (sync.value() == 0 && sync.stop() == 0));

    // Initialize sync counter
    sync.set_stop(count);

    if (callback != 0) {
      new_refresh_wait(id_refresh,callback);
    }
  } else {
    if (callback != 0) new_refresh_exit(*refresh);
  }
}

//----------------------------------------------------------------------
void Block::new_refresh_wait (int id_refresh, int callback)
{
  CHECK_ID(id_refresh);

  Refresh * refresh = cello::refresh(id_refresh);

  if (refresh->is_active()) {

    // make sure the callback parameter matches that in the refresh object

    ASSERT3("Block::new_refresh_wait()",
	   "Refresh[%d] mismatch between refresh callback %d and parameter callback %d",
	    id_refresh,callback, refresh->callback(),
	    (callback == refresh->callback()) );

    // make sure we aren't already in a "ready" state

    RefreshState & state = new_refresh_state_list_[id_refresh];

    ASSERT1("Block::new_refresh_wait()",
	   "Refresh[%d] not in 'active' state",
	    id_refresh,
	    (state == RefreshState::ACTIVE) );

    // tell refresh we're ready to start processing messages

    state = RefreshState::READY;

    // process any existing messages in the refresh message list

    Sync & sync = new_refresh_sync_list_[id_refresh];
    for (auto id_msg=0;
	 id_msg<new_refresh_msg_list_[id_refresh].size();
	 id_msg++) {

      MsgRefresh * msg = new_refresh_msg_list_[id_refresh][id_msg];

      // unpack message data into Block data
      msg->update(data());

      delete msg;
      sync.advance();
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

  Refresh * refresh    = cello::refresh(id_refresh);
  RefreshState & state = new_refresh_state_list_[id_refresh];
  Sync & sync          = new_refresh_sync_list_[id_refresh];

  ASSERT1("Block::new_refresh_check_done()",
	  "Refresh[%d] must not be in inactive state",
	  id_refresh,
	  (state != RefreshState::INACTIVE) );

  if (sync.stop()==0 || state == RefreshState::READY && (sync.is_done())) {

    // Make sure incoming message queue is empty

    ASSERT2("Block::new_refresh_wait()",
	   "Refresh %d message list has size %lu instead of 0",
	    id_refresh,new_refresh_msg_list_[id_refresh].size(),
	    (new_refresh_msg_list_[id_refresh].size() == 0));

    // reset sync counter
    sync.reset();
    sync.set_stop(0);

    // reset refresh state to inactive

    state = RefreshState::INACTIVE;

    // Call callback

    new_refresh_exit(*refresh);
  }
}

//----------------------------------------------------------------------

void Block::p_new_refresh_recv (MsgRefresh * msg)
{

  const int id_refresh = msg->id_refresh();
  CHECK_ID(id_refresh);

  RefreshState & state = new_refresh_state_list_[id_refresh];
  Sync & sync          = new_refresh_sync_list_[id_refresh];

  if (state == RefreshState::READY) {
    // unpack message data into Block data if ready
    msg->update(data());

    delete msg;

    sync.advance();

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

int Block::new_refresh_delete_particle_copies_ (Refresh * refresh){

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

  int count = 0;
  for (auto it_type=type_list.begin(); it_type != type_list.end(); it_type++){
    int it = *it_type;

    const int ia_c = particle.attribute_index(it,"is_local");
    const int cd   = particle.stride(it, ia_c);

    const int nb = particle.num_batches(it);

    int * is_local=0;

    for (int ib = 0; ib<nb; ib++){
      const int np = particle.num_particles(it,ib);
      is_local = (int *) particle.attribute_array(it, ia_c, ib);

      bool * mask = new bool[np];
      for( int ip=0; ip<np; ip++){
        mask[ip] = !(is_local[ip*cd]);
      }

      count += particle.delete_particles(it,ib,mask);

      delete [] mask;


    }
  }

  return count;
}

//----------------------------------------------------------------------


int Block::new_refresh_load_particle_copy_ (Refresh & refresh)
{
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

  int nl = particle_load_copy_
    (npa,particle_list,particle_array, index_list, &refresh);

  // Send particle data to neighbors

  new_particle_send_(refresh,nl,index_list,particle_list);

  delete [] index_list;


  return nl;
}

//----------------------------------------------------------------------

int Block::new_refresh_load_particle_faces_ (Refresh & refresh)
{
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

  //  TRACE_REFRESH("particle_load_faces()");

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

int Block::particle_load_copy_ (int npa,
				 ParticleData * particle_list[],
				 ParticleData * particle_array[],
				 Index index_list[],
				 Refresh *refresh)
{

  int nl = particle_create_array_neighbors_
    (refresh, particle_array,particle_list,index_list);

  // Scatter particles among particle_data array

  Particle particle (cello::particle_descr(),
		     data()->particle_data());

  std::vector<int> type_list;
  if (refresh->all_particles_copy()) {
    const int nt = particle.num_types();
    type_list.resize(nt);
    for (int i=0; i<nt; i++) type_list[i] = i;
  } else {
    type_list = refresh->particle_list_copy();
  }

  const bool copy = true;
  particle_scatter_neighbors_(npa,particle_array,type_list, particle,
                              copy);

  return nl;
}

//----------------------------------------------------------------------

int Block::particle_create_array_neighbors_
(Refresh * refresh,
 ParticleData * particle_array[],
 ParticleData * particle_list[],
 Index index_list[])
{
  //  TRACE_REFRESH("particle_create_array_neighbors()");

  const int rank = cello::rank();
  const int level = this->level();

  const int min_face_rank = refresh->min_face_rank();

  ItNeighbor it_neighbor =
    this->it_neighbor(min_face_rank,index_, neighbor_leaf,0,0);

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

  double dpx[nl],dpy[nl],dpz[nl];

  for (int i=0; i<nl; i++) {
    dpx[i]=0.0;
    dpy[i]=0.0;
    dpz[i]=0.0;
  }

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
 Particle particle,
 const bool copy  // default : false
 )
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
    const int ia_c  = particle.attribute_index(it, "is_local");

    // (...positions may use absolute coordinates (float) or
    // block-local coordinates (int))
    const bool is_float =
      (cello::type_is_float(particle.attribute_type(it,ia_x)));

    // (...stride may be != 1 if particle attributes are interleaved)
    const int d  = particle.stride(it,ia_x);

    //
    const int cd = particle.stride(it, ia_c);

    // ...for each batch of particles

    const int nb = particle.num_batches(it);

    int * is_local=0;

    for (int ib=0; ib<nb; ib++) {

      const int np = particle.num_particles(it,ib);

      if (np == 0) continue;

      // ...extract particle position arrays

      double * xa = new double [np];
      double * ya = new double [np];
      double * za = new double [np];

      particle.position(it,ib,xa,ya,za);

      is_local = (int *) particle.attribute_array(it, ia_c, ib);

      // ...initialize mask used for scatter and delete
      // ...and corresponding particle indices

      bool * mask = new bool[np];
      int  * index = new int[np];
      for (int ip=0; ip<np; ip++) {
// AJE - issue is here since in copying partilces arent on the grid boundaries
//    so mismatch between here and get_particle_bin_limits
//    which is called elsewhere in control_new_refresh
//    and now I realize this will currently only push copies of particles
//    fomr adjacent 1/2 of sibling grids. this is probably OK assuming
//    feedback region < 1/2 of grid size (probably?)
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
  if (copy){  // only copy particles that are not getting moved
    mask[ip] = in_block;
    is_local[ip*cd] = 0; // will not be local

    // hack - pretend partcles are closer to face
    ix = ix <= 1 ? 0 : 3;
    iy = iy <= 1 ? 0 : 3;
    iz = iz <= 1 ? 0 : 3;
    index[ip] = ix + 4 * (iy + 4 * iz);

  } else {    // only move particles that leave the block
    mask[ip] = ! in_block;
    is_local[ip*cd] = 1; // will be local once moved
  }
      }

      delete [] xa;
      delete [] ya;
      delete [] za;

      // ...scatter particles to particle array
      particle.scatter (it,ib, np, mask, index, npa, particle_array);
      // ... delete scattered particles if moved
      if (!copy) count += particle.delete_particles (it,ib,mask);


      delete [] mask;
      delete [] index;
    }
  }

  if (!copy) cello::simulation()->data_delete_particles(count);

}
