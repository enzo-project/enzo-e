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
  printf ("%s:%d TRACE_REFRESH %s\n",__FILE__,__LINE__,msg);
#else
#  define TRACE_REFRESH(msg) /* NOTHING */
#endif


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


  //  int ic3[3] = {0,0,0};

  if (refresh_type == refresh_coarse) {
    index_.child(index_.level(),ic3,ic3+1,ic3+2);
  }

  // ... copy field ghosts using FieldFace object.
  int n; char * array;
  bool lg3[3] = {false,false,false};
  std::vector<int> field_list = refresh()->field_list();

  FieldFace * field_face = load_field_face 
    (&n, &array, if3, ic3, lg3, refresh_type, field_list);

  // ... send the face data to the neighbor
  const int of3[3] = {-if3[0], -if3[1], -if3[2]};

  thisProxy[index_neighbor].p_refresh_store_field_face
    (n,array, refresh_type, of3, ic3);

  // ... delete the FieldFace created by load_face()
  delete field_face;
}

//----------------------------------------------------------------------

void Block::refresh_store_field_face_
(int n, char * buffer, int refresh_type, 
 int if3[3], int ic3[3])
{
  TRACE_REFRESH("refresh_store_field_face()");

  if (n > 0) {
    bool lg3[3] = {false,false,false};

    Refresh * refresh = this->refresh();
    std::vector<int> field_list = refresh->field_list();

    store_field_face (n,buffer, if3, ic3, lg3, refresh_type, field_list);
  }
}

//----------------------------------------------------------------------

void Block::refresh_load_particle_faces_ (Refresh *refresh)
{
  const int min_face_rank = refresh->min_face_rank();
  const int neighbor_type = refresh->neighbor_type();

  const int rank = this->rank();
  int np = rank == 1 ? 4 : (rank == 2 ? 4*4 : 4*4*4);

  // 1. CREATE PARTICLE_DATA ARRAY

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
  // |       ║   :   ║    
  // |   2   + - + - +   3
  // |       ║   :   ║    
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

  const int level = this->level();

  ParticleData * particle_data[np];

  // 2. SCATTER PARTICLES AMONG PARTICLE_DATA ARRAY

  ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);

  int count = 0;
  while (it_neighbor.next()) {

    Index index_neighbor = it_neighbor.index();

    const int level_face = it_neighbor.face_level();

    int if3[3] = {0,0,0} ,ic3[3] = {0,0,0};

    it_neighbor.face(if3);

    const int refresh_type = 
      (level_face == level - 1) ? refresh_coarse :
      (level_face == level)     ? refresh_same :
      (level_face == level + 1) ? refresh_fine : refresh_unknown;

    if (refresh_type==refresh_coarse) {
      index_.child(index_.level(),ic3,ic3+1,ic3+2);
    } else if (refresh_type==refresh_fine) {
      it_neighbor.child(ic3);
    }

    int lower[3] = {0,0,0};
    int upper[3] = {1,1,1};
    refresh->index_limits (rank,refresh_type,if3,ic3,lower,upper);

    for (int ix=lower[0]; ix<upper[0]; ix++) {
      for (int iy=lower[1]; iy<upper[1]; iy++) {
	for (int iz=lower[2]; iz<upper[2]; iz++) {
	  count ++;
	}
      }
    }
    printf ("refresh %s face %d %d %d child %d %d %d lower %d %d %d upper %d %d %d\n",
	    refresh_type == refresh_coarse ? "coarse" :
	    ( refresh_type == refresh_same ? "same" : "fine"),
	    if3[0],if3[1],if3[2],
	    ic3[0],ic3[1],ic3[2],
	    lower[0],lower[1],lower[2],
	    upper[0],upper[1],upper[2]);

  }

  printf ("refresh count = %d\n",count);
  // if (neighbor_type == neighbor_leaf) {


  // } else if (neighbor_type == neighbor_level) {

  //   ERROR ("Block::refresh_load_particle_faces_()",
  // 	   "neighbor_level not supported (designed for MG)");
  // }
}

//----------------------------------------------------------------------

void Block::refresh_load_particle_face_
( int refresh_type,
  Index index_neighbor,
  int if3[3],
  int ic3[3])

{
  TRACE_REFRESH("refresh_load_particle_face()");

  std::vector<int> particle_list = refresh()->particle_list();

  Particle particle (simulation()->particle_descr(),
		     data()->particle_data());

  const int nt = particle_list.size();

  // number of calls to p_refresh_store_particle_face

  for (int it=0; it<nt; it++) {

    const bool interleaved = particle.interleaved(it);
    const int nb = particle.num_batches(it);

    for (int ib=0; ib<nb; ib++) {

      const int np = particle.num_particles(it,ib);
      const int mb = particle.batch_size();
      const int mp = particle.particle_bytes(it);

      // if not interleaved, batch size is always full

      int n = mp*(interleaved ? np : mb);
      char * array = particle.attribute_array(it,0,ib);

      thisProxy[index_neighbor].p_refresh_store_particle_face
	(n,array, it, ib, if3, ic3 );
    }
  }
}

//----------------------------------------------------------------------

void Block::refresh_store_particle_face_
(int n, char * array, 
 int it, int ib,
 int if3[3], 
 int ic3[3])
{
  TRACE_REFRESH("refresh_store_particle_face");

  Particle particle (simulation()->particle_descr(),
		     data()->particle_data());

  // determine number of particles

  int mp = particle.particle_bytes(it);
  const int mb = particle.batch_size();

  const bool interleaved = particle.interleaved(it);
  int np = (interleaved) ? n/mp : mb;

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
      for (int iy=0; iy<ny; iy++) {
	a_dst [iy + mp*jp] = array [iy + mp*ip];
      }
    }
    jp = (jp+1) % mb;
    if (jp==0) jb++;
  }

}
