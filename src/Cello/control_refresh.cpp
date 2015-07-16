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

//----------------------------------------------------------------------

void Block::refresh_begin_() 
{

  Refresh * refresh = this->refresh();

  check_leaf_();

  check_delete_();

  refresh->sync_load().reset();

  simulation()->set_phase(phase_refresh);

  // Refresh if Refresh object exists and have data

  int if3[3];
  int ic3[3];

  int count = 0;

  if ( refresh && refresh->active() ) {

    const int min_face_rank = refresh->min_face_rank();

    const int level = this->level();

    if (refresh->neighbor_type() == neighbor_leaf) {
      ItNeighbor it_neighbor = this->it_neighbor(min_face_rank,index_);
      int if3[3];
      int ic3[3];
      while (it_neighbor.next()) {
	Index index_neighbor = it_neighbor.index();
	it_neighbor.face(if3);
	it_neighbor.child(ic3);
	int level_face = it_neighbor.face_level();
	++count;

	if (level_face == level) {
	  refresh_load_face_ (refresh_same,index_neighbor,if3,ic3);
	} else if (level_face == level + 1) {
	  refresh_load_face_ (refresh_fine,index_neighbor,if3,ic3);
	} else if (level_face == level - 1) {
	  refresh_load_face_ (refresh_coarse,index_neighbor,if3,ic3);
	}

      }
    } else if (refresh->neighbor_type() == neighbor_level) {
      ItFace it_face = this->it_face(min_face_rank,index_);
      int if3[3];
      int ic3[3] = {0,0,0}; // not accessed since same level
      while (it_face.next()) {
	Index index_face = it_face.index();
	it_face.face(if3);
	++count;

	refresh_load_face_ (refresh_same,index_face,if3,ic3);

      }
    }
  }

  // call with self to set counter
  refresh_load_face_(refresh_same,index(),if3,ic3,count + 1);

}

//----------------------------------------------------------------------

void Block::refresh_load_face_
( int type_refresh,
  Index index_neighbor,
  int iface[3],
  int ichild[3],
  int count)
{
  Refresh * refresh = this->refresh();

  if (count != 0) {

     refresh->sync_load().set_stop(count);

  } else {

    FieldFace * field_face;

    int type_op_array;

    int level = index_.level();
    switch (type_refresh) {
    case refresh_coarse:

      index_.child(level,ichild,ichild+1,ichild+2);

      type_op_array = op_array_restrict;

      break;

    case refresh_same:

      type_op_array = op_array_copy;

      break;

    case refresh_fine:

      type_op_array = op_array_prolong;

      break;

    default:
      type_op_array = op_array_unknown;
      ERROR1("Block::refresh_face()",
	     "Unknown type_refresh %d",
	     type_refresh);

    }

    int n; 
    char * array;
    bool lghost[3] = {false,false,false};

    std::vector<int> field_list = refresh->field_list();

    field_face = load_face (&n, &array,
			    iface, ichild, lghost,
			    type_op_array,
			    field_list);

    int jface[3] = {-iface[0], -iface[1], -iface[2]};

    thisProxy[index_neighbor].p_refresh_store_face
      (n,array, type_refresh, jface, ichild);

    delete field_face;
  }

  if (refresh->sync_load().next()) {
    
    control_sync (CkIndex_Block::p_refresh_exit(),refresh_.sync_type(),2);
  }
}

//----------------------------------------------------------------------

void Block::refresh_store_face_
(int n, char * buffer, int type_refresh,
 int iface[3], int ichild[3],int count
)
{
  Refresh * refresh = this->refresh();

  if (count==0) {
    bool lghost[3] = {false,false,false};

    std::vector<int> field_list = refresh->field_list();

    int op_array;
    switch (type_refresh) {
    case refresh_coarse:  op_array = op_array_restrict;  break;
    case refresh_same:    op_array = op_array_copy;      break;
    case refresh_fine:    op_array = op_array_prolong;   break;
    default:
      op_array = op_array_unknown;
      ERROR1("Block::refresh_store_face_()",
	     "Unknown type_refresh %d",  type_refresh);
      break;
    }

    store_face_(n,buffer,
		iface, ichild, lghost,
		op_array,
		field_list);
  }
}

//----------------------------------------------------------------------
