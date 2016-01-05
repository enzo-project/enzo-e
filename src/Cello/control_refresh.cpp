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

/* #define DEBUG_REFRESH */

#ifdef DEBUG_REFRESH
#  define TRACE_REFRESH(msg) \
  printf ("%s:%d %p TRACE_REFRESH %s\n",__FILE__,__LINE__,this,msg);	\
  fflush(stdout);
#else
#  define TRACE_REFRESH(msg) /* NOTHING */
#endif

/* #define NEW_REFRESH */

//----------------------------------------------------------------------

void Block::refresh_begin_() 
{
  TRACE_REFRESH("refresh_begin_()");

  Refresh * refresh = this->refresh();

  check_leaf_();

  check_delete_();

  simulation()->set_phase(phase_refresh);

  // Refresh if Refresh object exists and have data

  if ( refresh && refresh->active() ) {

    refresh_load_field_faces_   (refresh);
    refresh_load_particle_faces_(refresh);

  }

  control_sync (CkIndex_Block::p_refresh_exit(),refresh_.sync_type(),2);
}

//----------------------------------------------------------------------

void Block::refresh_load_field_faces_ (Refresh *refresh)
{
  const int min_face_rank = refresh->min_face_rank();
  const int neighbor_type = refresh->neighbor_type();

  if (neighbor_type == neighbor_leaf) {

    // Loop over neighbor leaf Blocks (not necessarily same level)

    ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);

    while (it_neighbor.next()) {

      Index index_neighbor = it_neighbor.index();

      int if3[3],ic3[3];
      it_neighbor.face (if3);
      it_neighbor.child(ic3);

      const int level = this->level();
      const int level_face = it_neighbor.face_level();

      const int refresh_type = 
	(level_face == level - 1) ? refresh_coarse :
	(level_face == level)     ? refresh_same :
	(level_face == level + 1) ? refresh_fine : refresh_unknown;

      refresh_load_field_face_ (refresh_type,index_neighbor,if3,ic3);
    }

  } else if (neighbor_type == neighbor_level) {

    // Loop over neighbor Blocks in same level (not necessarily leaves)

    ItFace it_face = this->it_face(min_face_rank,index_);

    while (it_face.next()) {

      Index index_face = it_face.index();

      int if3[3],ic3[3];
      it_face.face(if3);

      refresh_load_field_face_ (refresh_same,index_face,if3,ic3);

    }
  }
}

//----------------------------------------------------------------------

void Block::refresh_load_field_face_
( int refresh_type,
  Index index_neighbor,
  int if3[3],
  int ic3[3])

{
  TRACE_REFRESH("refresh_load_field_face()");

  // REFRESH FIELDS

  // ... coarse neighbor requires child index of self in parent

  if (refresh_type == refresh_coarse) {
    index_.child(index_.level(),ic3,ic3+1,ic3+2);
  }

  // ... copy field ghosts to array using FieldFace object
  int n; char * array;
  bool lg3[3] = {false,false,false};
  std::vector<int> field_list = refresh()->field_list();

  FieldFace * field_face = load_field_face 
    (&n, &array, if3, ic3, lg3, refresh_type, field_list);


#ifdef NEW_REFRESH

  DataMsg * data_msg = new DataMsg;
  data_msg -> set_field_face (field_face);
  thisProxy[index_neighbor].p_refresh_store (data_msg);

#else

  // ... send the face data to the neighbor
  const int of3[3] = {-if3[0], -if3[1], -if3[2]};

  thisProxy[index_neighbor].p_refresh_store_field_face
    (n,array, refresh_type, of3, ic3);

#endif

  // ... delete the FieldFace created by load_face()
  delete field_face;
  delete [] array;
}


//----------------------------------------------------------------------


void Block::p_refresh_store (DataMsg * msg)
{
#ifdef NEW_REFRESH
 msg->update(data());
#endif
}


//----------------------------------------------------------------------


void Block::refresh_store_field_face_
(int n, char * array, int refresh_type, 
 int if3[3], int ic3[3])
{
#ifndef NEW_REFRESH
  TRACE_REFRESH("refresh_store_field_face()");

  if (n > 0) {

    // copy array to FieldData

    bool lg3[3] = {false,false,false};

    Refresh * refresh = this->refresh();
    std::vector<int> field_list = refresh->field_list();

    store_field_face (n,array, if3, ic3, lg3, refresh_type, field_list);
  }
#endif
}


//----------------------------------------------------------------------

void Block::refresh_load_particle_faces_ (Refresh *refresh)
{
  // Array elements correspond to child-sized blocks to
  // the left, inside, and right of the main Block.  Particles
  // are assumed to be (well) within this area.
  //
  //     +---+---+---+---+
  //     | 03| 13| 23| 33|
  //     +---+===+===+---+
  //     | 02║   :   ║ 32|
  //     +---+ - + - +---+
  //     | 01║   :   ║ 31|
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
  // |       ║       ║    
  // |   2   +       +   3
  // |       ║       ║    
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
  //     | 2 ║   :   ║ 3 |
  //     +---+ - + - +---+
  //     | 2 ║   :   ║ 3 |
  //     +---+=======+---+
  //     | 0 | 1 | 1 | 1 |
  //     +---+---+---+---+

  TRACE_REFRESH("refresh_load_particle_faces()");

  const int min_face_rank = refresh->min_face_rank();

  const int level = this->level();
  const int rank = this->rank();

  // 1. CREATE PARTICLE_DATA ARRAY

  // ... size of the ParticleData array
  const int npa = ( (rank == 1) ? 4 : ( (rank == 2) ? 4*4 : 4*4*4));

  ParticleData * particle_array[npa];
  ParticleData * particle_list[npa];
  Index index_list[npa];
  
  for (int i=0; i<npa; i++) {
    particle_list[i] = NULL;
    particle_array[i] = NULL;
  }

  // 2. INITIALIZE THE PARTICLE_DATA ARRAY

  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);

  // ... prepare for adjusting position based on periodicity

  //     ... domain extents

  double dxm,dym,dzm;
  double dxp,dyp,dzp;

  simulation()->hierarchy()->lower(&dxm,&dym,&dzm);
  simulation()->hierarchy()->upper(&dxp,&dyp,&dzp);

  //     ... periodicity
  bool p32[3][2];
  periodicity(p32);

  //     ... boundary
  bool b32[3][2];
  is_on_boundary (b32);

  //     ... position adjust for each particle array
  double dpx[npa],dpy[npa],dpz[npa];
  for (int i=0; i<npa; i++) { dpx[i]=dpy[i]=dpz[i]=0.0; }

  int count = 0;
  int nl = 0;
  while (it_neighbor.next()) {

    const int level_face = it_neighbor.face_level();

    int if3[3] = {0,0,0} ,ic3[3] = {0,0,0};

    it_neighbor.face(if3);

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

    int lower[3] = {0,0,0};
    int upper[3] = {1,1,1};
    refresh->index_limits (rank,refresh_type,if3,ic3,lower,upper);

    ParticleData * pd = new ParticleData;
    ParticleDescr * particle_descr = simulation() -> particle_descr();
    pd->allocate(particle_descr);

    particle_list[nl] = pd;
    index_list[nl] = it_neighbor.index();

    // ... adjust for periodic boundary conditions
    if (rank >= 1) {
      if (lower[0]==0 && b32[0][0] && p32[0][0]) dpx[nl] = +(dxp - dxm);
      if (upper[0]==4 && b32[0][1] && p32[0][1]) dpx[nl] = -(dxp - dxm);
    }
    if (rank >= 2) {
      if (lower[1]==0 && b32[1][0] && p32[1][0]) dpy[nl] = +(dyp - dym);
      if (upper[1]==4 && b32[1][1] && p32[1][1]) dpy[nl] = -(dyp - dym);
    }
    if (rank >= 3) {
      if (lower[2]==0 && b32[2][0] && p32[2][0]) dpz[nl] = +(dzp - dzm);
      if (upper[2]==4 && b32[2][1] && p32[2][1]) dpz[nl] = -(dzp - dzm);
    }

    ++ nl;

    for (int iz=lower[2]; iz<upper[2]; iz++) {
      for (int iy=lower[1]; iy<upper[1]; iy++) {
	for (int ix=lower[0]; ix<upper[0]; ix++) {
	  int i=ix + 4*(iy + 4*iz);
	  particle_array[i] = pd;
	  count ++;
	}
      }
    }
  }

  // 3. SCATTER PARTICLES AMONG PARTICLE_DATA ARRAY

  std::vector<int> type_list = refresh->particle_list();


  Particle particle (simulation() -> particle_descr(),
		     data()       -> particle_data());
  //     ... get Block bounds 
  double xm,ym,zm;
  lower(&xm,&ym,&zm);
  double xp,yp,zp;
  upper(&xp,&yp,&zp);

  // find block center (x0,y0,z0) and width (xl,yl,zl)
  const double x0 = 0.5*(xm+xp);
  const double y0 = 0.5*(ym+yp);
  const double z0 = 0.5*(zm+zp);
  const double xl = xp-xm;
  const double yl = yp-ym;
  const double zl = zp-zm;

  // ...for each particle type to be moved
  std::vector<int>::iterator it_type;
  for (it_type=type_list.begin(); it_type!=type_list.end(); it_type++) {

    int it = *it_type;

    const int ia_x = particle.attribute_position(it,0);
    const int ia_id = particle.attribute_index(it,"id");
    const bool is_float = (cello::type_is_float(particle.attribute_type(it,ia_x)));

    const int d = particle.stride(it,ia_x);
    const int di = particle.stride(it,ia_id);

    const int nb = particle.num_batches(it);

    // ...for each batch of particles
    for (int ib=0; ib<nb; ib++) {

      const int np = particle.num_particles(it,ib);

      // ...extract particle positions

      double xa[np],ya[np],za[np];
      particle.position(it,ib,xa,ya,za);

      int64_t * ida = (int64_t*)particle.attribute_array (it,ia_id,ib);
      // ...initialize masks (scatter mask and delete mask) and index array
      bool mask_scatter[np],mask_delete[np];
      int index[np];
      for (int ip=0; ip<np; ip++) {

	int64_t id = ida[ip*di];
	double x = is_float ? 2.0*(xa[ip*d]-x0)/xl : xa[ip*d];
	double y = is_float ? 2.0*(ya[ip*d]-y0)/yl : ya[ip*d];
	double z = is_float ? 2.0*(za[ip*d]-z0)/zl : za[ip*d];

	int ix = (rank >= 1) ? (x + 2) : 0;
	int iy = (rank >= 2) ? (y + 2) : 0;
	int iz = (rank >= 3) ? (z + 2) : 0;
	int i = ix + 4*(iy + 4*iz);

	if (! (0 <= ix && ix < 4) ) {
	  printf ("id %ld\n",id);
	  printf ("ix iy iz %d %d %d\n",ix,iy,iz);
	  printf ("x y z %f %f %f\n",x,y,z);
	  printf ("xa ya za %f %f %f\n",xa[ip*d],ya[ip*d],za[ip*d]);
	  printf ("xm ym zm %f %f %f\n",xm,ym,zm);
	  printf ("xp yp zp %f %f %f\n",xp,yp,zp);
	  fflush(stdout);
	  ASSERT1 ("Block::refresh_load_particle_faces_",
		   "particle data index ix = %d out of bounds",
		   ix,(0 <= ix && ix < 4));
	}
	ASSERT1 ("Block::refresh_load_particle_faces_",
		 "particle data index iy = %d out of bounds",
		 iy,(0 <= iy && iy < 4));
	fflush(stdout);
	ASSERT1 ("Block::refresh_load_particle_faces_",
		 "particle data index iz = %d out of bounds",
		 iz,(0 <= iz && iz < 4));
	fflush(stdout);
	index[ip] = i;
	bool in_block = true;
	in_block = in_block && (!(rank >= 1) || (1 <= ix && ix <= 2));
	in_block = in_block && (!(rank >= 2) || (1 <= iy && iy <= 2));
	in_block = in_block && (!(rank >= 3) || (1 <= iz && iz <= 2));
	mask_scatter[ip] = ! in_block;
	mask_delete[ip] = mask_scatter[ip];
      }

      // ...scatter particles to particle array
      particle.scatter (it,ib, np, mask_scatter, index, npa, particle_array);
      // ... delete scattered particles
      particle.delete_particles (it,ib,mask_delete);
    }
  }

  // 3.5 Update positions x += dpx[i], etc.

  // ... loop over index list not array to prevent multiple updates

  ParticleDescr * p_descr = simulation() -> particle_descr();

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

	  const int np = particle_neighbor.num_particles(it,ib);
	  double xa[np],ya[np],za[np];

	  particle_neighbor.position(it,ib,xa,ya,za);

	  particle_neighbor.position_update (it,ib,dpx[il],dpy[il],dpz[il]);

	}
      }
    }
  }

  // 4. SEND PARTICLE DATA TO NEIGHBORS

  for (int il=0; il<nl; il++) {

    Index index_neighbor  =    index_list[il];
    ParticleData * p_data = particle_list[il];
    Particle particle_send (p_descr,p_data);
    
    const int nt = particle_send.num_types();

    for (int it=0; it<nt; it++) {

      const bool interleaved = particle_send.interleaved(it);
      const int nb = particle_send.num_batches(it);

      for (int ib=0; ib<nb; ib++) {

	const int np = particle_send.num_particles(it,ib);
	const int mb = particle_send.batch_size();
	int mp = particle_send.particle_bytes(it);

	// if not interleaved, batch size is always full

	int n = mp*(interleaved ? np : mb);
	char * a = particle_send.attribute_array(it,0,ib);

	thisProxy[index_neighbor].p_refresh_store_particle_face (n,np,a,it);
      }
    }
    
  }


  // delete neighbor particle objects
  for (int il=0; il<nl; il++) {
    ParticleData * p_data = particle_list[il];
    delete p_data;
  }
}

//----------------------------------------------------------------------

void Block::refresh_store_particle_face_
  (int n, int np, char * a, int it)
{
  TRACE_REFRESH("refresh_store_particle_face");

  Particle particle (simulation() -> particle_descr(),
		     data()       -> particle_data());

  // determine number of particles

  int mp = particle.particle_bytes(it);
  const int mb = particle.batch_size();

  const bool interleaved = particle.interleaved(it);

  int j = particle.insert_particles(it,np);

  int jb,jp;
  particle.index(j,&jb,&jp);

  // copy particles from array

  const int na = particle.num_attributes(it);

  for (int ip=0; ip<np; ip++) {
    for (int ia=0; ia<na; ia++) {
      if (!interleaved) 
	mp = particle.attribute_bytes(it,ia);
      int ny = particle.attribute_bytes(it,ia);
      char * a_dst = particle.attribute_array(it,ia,jb);
      int increment = particle.attribute_array(it,ia,jb) 
	-             particle.attribute_array(it,0,jb);
      char * a_src = a + increment;
      for (int iy=0; iy<ny; iy++) {
        a_dst [iy + mp*jp] = a_src [iy + mp*ip];
      }
    }

    jp = (jp+1) % mb;
    if (jp==0) jb++;
  }
}
