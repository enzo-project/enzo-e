// See LICENSE_CELLO file for license and copyright information

/// @file     control_refresh.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-26
/// @brief    Charm-related functions associated with refreshing ghost zones

#include "simulation.hpp"
#include "mesh.hpp"
#include "comm.hpp"

#include "charm_simulation.hpp"
#include "charm_mesh.hpp"

static char buffer[256];

//----------------------------------------------------------------------

void CommBlock::refresh_begin_() 
{

  Simulation * simulation = proxy_simulation.ckLocalBranch();

  const Config * config = simulation->config();

  if (is_leaf()) {

    int rank = simulation->dimension();

    int n3[3];
    size_forest(&n3[0],&n3[1],&n3[2]);

    int refresh_rank = config->field_refresh_rank;
    loop_refresh_.set_stop(0);

    ItFace it_face(rank,refresh_rank);
    int if3[3];

    const int level = this->level();

    while (it_face.next(if3)) {


      Index index_neighbor = index_.index_neighbor(if3,n3);

      TRACE0;

      if (face_level(if3) == level-1) {       // COARSE

	int ic3[3];
	index_.child(level,ic3+0,ic3+1,ic3+2);
	int ip3[3];
	parent_face_(ip3,if3,ic3);

	refresh_face
	  (refresh_coarse,index_neighbor.index_parent(),ip3,ic3);

      } else if (face_level(if3) == level) {    // SAME

	int ic3[3] = {0,0,0};
	refresh_face (refresh_same,index_neighbor,if3,ic3);

      } else if (face_level(if3) == level+1) {  // FINE
	    
	int ic3m[3];
	int ic3p[3];
	loop_limits_nibling_(ic3m,ic3p,if3);
	int ic3[3];
	for (ic3[0]=ic3m[0]; ic3[0]<=ic3p[0]; ic3[0]++) {
	  for (ic3[1]=ic3m[1]; ic3[1]<=ic3p[1]; ic3[1]++) {
	    for (ic3[2]=ic3m[2]; ic3[2]<=ic3p[2]; ic3[2]++) {

	      int jc3[3];
	      facing_child_ (jc3,ic3,if3);

	      Index index_nibling = 
		index_neighbor.index_child(jc3[0],jc3[1],jc3[2]);
		  
	      refresh_face (refresh_fine,index_nibling, if3,ic3);
	    }
	  }
	}
      } else {
	sprintf (buffer,"REFRESH ERROR face (%d %d %d) level %d face_level %d phase %d",
		 if3[0],if3[1],if3[2],level,face_level(if3),next_phase_);
	index_.print(buffer,-1,2,false,simulation);
      
	ERROR("CommBlock::refresh_begin_()",
	      "Refresh error");
      }

    }
  }

  control_sync (sync_refresh_exit);

}

//----------------------------------------------------------------------

void CommBlock::refresh_face
( int type_refresh,
  Index index_neighbor,
  int iface[3],
  int ichild[3] )
{

  int n; 
  char * array;

  bool lghost[3] = {false,false,false};

  FieldFace * field_face;

  int type_op_array;

  if (type_refresh == refresh_coarse) {

    int level = index_.level();
    index_.child(level,ichild,ichild+1,ichild+2);

    type_op_array = op_array_restrict;

  } else if (type_refresh == refresh_same) {

    type_op_array = op_array_copy;


  } else if (type_refresh == refresh_fine) {

    type_op_array = op_array_prolong;

  }

  field_face = load_face_ (&n, &array,
			   iface, ichild, lghost,
			   type_op_array);

  int jface[3] = {-iface[0], -iface[1], -iface[2]};

#ifdef CELLO_DEBUG
  // std::string bit_str = index_neighbor.bit_string(-1,2);
  // sprintf (buffer," -> %s",bit_str.c_str());
  // index_.print(buffer,-1,2,false,simulation());
#endif  

  // --------------------------------------------------
  // ENTRY: #2 CommBlock::refresh_face()-> CommBlock::x_refresh_face()
  // ENTRY: neighbor
  // --------------------------------------------------
  thisProxy[index_neighbor].x_refresh_face
    (n,array, type_refresh, jface, ichild);
  // --------------------------------------------------

  loop_refresh_.add_stop();
  delete field_face;
}

//----------------------------------------------------------------------

void CommBlock::refresh_face_ (int n, char * buffer, int type_refresh,
			       int iface[3], int ichild[3])
{
  if (type_refresh == refresh_coarse) { // coarse

    bool lghost[3] = {false};

    store_face_(n,buffer,
		iface, ichild, lghost,
		op_array_restrict);

  } else if (type_refresh == refresh_same) { // same

    if ( n != 0) {
      bool lghost[3] = {false,false,false};
      store_face_(n,buffer,
		  iface,ichild,lghost,
		  op_array_copy);
    }

  } else if (type_refresh == refresh_fine) {

    bool lghost[3] = {false};
    store_face_(n,buffer,
		iface, ichild, lghost,
		op_array_prolong);
  }
}

//----------------------------------------------------------------------

void CommBlock::x_refresh_child 
(
 int    n, 
 char * buffer, 
 int    ic3[3]
 )
{
  int  iface[3]  = {0,0,0};
  bool lghost[3] = {true,true,true};
  store_face_(n,buffer, iface, ic3, lghost, op_array_restrict);
}

//----------------------------------------------------------------------
