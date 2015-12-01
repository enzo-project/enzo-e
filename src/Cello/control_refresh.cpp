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

  refresh->sync_load().reset();

  simulation()->set_phase(phase_refresh);

  // Refresh if Refresh object exists and have data

  int if3[3] = {0,0,0};
  int ic3[3] = {0,0,0};

  int count = 0;

  if ( refresh && refresh->active() ) {

    const int min_face_rank = refresh->min_face_rank();
    const int neighbor_type = refresh->neighbor_type();

    if (neighbor_type == neighbor_leaf) {

      ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);

      while (it_neighbor.next()) {

	Index index_neighbor = it_neighbor.index();

	// (shadowing input parameters!)
	int if3[3];
	int ic3[3];

	it_neighbor.face (if3);
	it_neighbor.child(ic3);

	const int level = this->level();
	const int level_face = it_neighbor.face_level();
	const int refresh_type = 
	  (level_face == level - 1) ? refresh_coarse :
	  (level_face == level)     ? refresh_same :
	  (level_face == level + 1) ? refresh_fine : refresh_unknown;

	count += refresh_load_face_ (refresh_type,index_neighbor,if3,ic3);
      }

    } else if (neighbor_type == neighbor_level) {

      ItFace it_face = this->it_face(min_face_rank,index_);
      while (it_face.next()) {
	Index index_face = it_face.index();
	it_face.face(if3);

	count += refresh_load_face_ (refresh_same,index_face,if3,ic3);

      }
    }
  }

  // call with self to set counter
  refresh_load_face_(refresh_same,index(),if3,ic3,count + 1);

}

//----------------------------------------------------------------------

int Block::refresh_load_face_
( int refresh_type,
  Index index_neighbor,
  int if3[3],
  int ic3[3],
  int count)
{
  TRACE_REFRESH("refresh_load_face()");

  Sync & sync_load = this->refresh()->sync_load();

  int num_msg = 0;

  if (count != 0) {

     sync_load.inc_stop(count);

  } else {

    num_msg += refresh_load_field_face_   
      (refresh_type,index_neighbor,if3,ic3);

    num_msg += refresh_load_particle_face_
      (refresh_type,index_neighbor,if3,ic3);

  }

  // If all faces updated, exit refresh

  if (sync_load.next()) {
    control_sync (CkIndex_Block::p_refresh_exit(),refresh_.sync_type(),2);
  }

  return num_msg;
}

//----------------------------------------------------------------------

void Block::refresh_store_face_
(int n, char * buffer, int refresh_type,
 int if3[3], int ic3[3],int count
 )
{
  TRACE_REFRESH("refresh_store_face()");
  if (count==0) {
    refresh_store_field_face_   (0,buffer,refresh_type,if3,ic3);
    refresh_store_particle_face_(0,buffer,0,0,if3,ic3);
  }
}

//----------------------------------------------------------------------

int Block::refresh_load_field_face_
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

  return 1;
}

//----------------------------------------------------------------------

void Block::refresh_store_field_face_
(int n, char * buffer, int refresh_type, 
 int if3[3], int ic3[3],int count )
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

int Block::refresh_load_particle_face_
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

  int count = 0;

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

      ++ count;
	
				     
    }
  }
  return count;
}

//----------------------------------------------------------------------

void Block::refresh_store_particle_face_
(int n, char * array, 
 int it, int ib,
 int if3[3], 
 int ic3[3], 
 int count )
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
