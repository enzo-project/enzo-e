// See LICENSE_CELLO file for license and copyright information

/// @file     control_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with refreshing ghost zones
/// @ingroup  Control

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"
#include "control.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

// static char buffer[256];

//----------------------------------------------------------------------


void CommBlock::refresh_begin_() 
{
  Refresh * refresh = simulation()->problem()->refresh(index_refresh_);

  // Check consistency between is_leaf() and std::vector children_
  if ((  is_leaf() && children_.size() != 0) ||
      (! is_leaf() && children_.size() == 0)) {
    const char * logic_str[2] = {"false","true"};
    WARNING4("CommBlock::refresh_begin_()",
	     "%s: is_leaf() == %s && children_.size() == %d"
	     "setting is_leaf_ <== %s",
	     name_.c_str(), logic_str[is_leaf()?1:0],
	     children_.size(),logic_str[is_leaf()?0:1]);
    is_leaf_ = ! is_leaf();
  }

  // Check if CommBlock should have been deleted
  if (delete_) {
    WARNING1("refresh_begin_()",
	     "%s: refresh called on deleted CommBlock",
	     name_.c_str());
    return;
  }

  simulation()->set_phase(phase_refresh);

  // Refresh if Refresh object exists and have data

  if (refresh && is_leaf()) {

    // Level of this CommBlock

    const int level = this->level();

    // Get forest size required for non-periodic boundaries

    int n3[3];
    size_forest(&n3[0],&n3[1],&n3[2]);

    // Get the face iterator for specified field face rank 
    // (0 corner, 1 edge, 2 facet)

    ItFace it_face = this->it_face( refresh->field_face_rank() );

    // Loop over specified faces

    int if3[3];
    while (it_face.next(if3)) {

      Index index_neighbor = index_.index_neighbor(if3,n3);

      TRACE0;

      if (face_level(if3) == level-1) {

	// If COARSE LEVEL neighbor, refresh neighbor's parent
	
	int ic3[3];
	index_.child(level,ic3+0,ic3+1,ic3+2);
	int ip3[3];
	parent_face_(ip3,if3,ic3);

	refresh_load_face_
	  (refresh_coarse,index_neighbor.index_parent(),ip3,ic3);

      } else if (face_level(if3) == level) {

	// if SAME LEVEL neighbor, refresh neighbor

	int ic3[3] = {0,0,0};
	refresh_load_face_ (refresh_same,index_neighbor,if3,ic3);

      } else if (face_level(if3) == level+1) {
	    
	// FINE LEVEL neighbor, refresh adjacent neighbor's children

	// get lower and upper limits of MY children along face if3
	int ic3m[3],ic3p[3];
	loop_limits_nibling_(ic3m,ic3p,if3);

	int ic3[3];
	for (ic3[0]=ic3m[0]; ic3[0]<=ic3p[0]; ic3[0]++) {
	  for (ic3[1]=ic3m[1]; ic3[1]<=ic3p[1]; ic3[1]++) {
	    for (ic3[2]=ic3m[2]; ic3[2]<=ic3p[2]; ic3[2]++) {

	      // get corresponding index of child in neighbor
	      int jc3[3];
	      facing_child_ (jc3,ic3,if3);

	      Index index_nibling = 
		index_neighbor.index_child(jc3[0],jc3[1],jc3[2]);
		  
	      refresh_load_face_ (refresh_fine,index_nibling, if3,ic3);
	    }
	  }
	}
      } else {

	index_.print("ERROR Refresh",-1,2,false,simulation());
	debug_faces_("refresh");
	ERROR7("CommBlock::refresh_begin_()",
	       "%s: REFRESH ERROR: "
	       "face (%d %d %d) "
	       "level %d "
	       "face_level %d "
	       "is_leaf %d",
	       name_.c_str(),
	       if3[0],if3[1],if3[2],
	       level,
	       face_level(if3),
	       is_leaf());
      }

    }
  }

  // WARNING("CommBlock::refresh_begin_",
  // 	  "refresh synchronization changed from contribute to neighbor"
  //      "(see bug #44)");

  control_next (phase_refresh_exit,"neighbor");
  //  control_next (phase_refresh_exit,"contribute");

}

//----------------------------------------------------------------------

void CommBlock::refresh_load_face_
( int type_refresh,
  Index index_neighbor,
  int iface[3],
  int ichild[3] )
{
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
    ERROR1("CommBlock::refresh_face()",
	   "Unknown type_refresh %d",
	   type_refresh);

  }

  int n; 
  char * array;
  bool lghost[3] = {false,false,false};

  std::vector<int> field_list;
  field_face = load_face_ (&n, &array,
			   iface, ichild, lghost,
			   type_op_array,
			   field_list);

  int jface[3] = {-iface[0], -iface[1], -iface[2]};

  thisProxy[index_neighbor].x_refresh_send_face
    (n,array, type_refresh, jface, ichild);

  delete field_face;
}

//----------------------------------------------------------------------

void CommBlock::refresh_store_face_
(int n, char * buffer, int type_refresh,
 int iface[3], int ichild[3])
{
  bool lghost[3] = {false,false,false};

  std::vector<int> field_list = refresh()->field_list();

  int op_array;
  switch (type_refresh) {
  case refresh_coarse:  op_array = op_array_restrict;  break;
  case refresh_same:    op_array = op_array_copy;      break;
  case refresh_fine:    op_array = op_array_prolong;   break;
  default:
    ERROR1("CommBlock::refresh_store_face_()",
	   "Unknown type_refresh %d",  type_refresh);
    break;
  }

  store_face_(n,buffer,
	      iface, ichild, lghost,
	      op_array,
	      field_list);
}

//----------------------------------------------------------------------
