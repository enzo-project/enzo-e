//======================================================================
//
//        File: hypre.cpp
//
//     Summary: Interface routines to HYPRE
//
// Description:Interface routines to HYPRE
//
//      Author: James Bordner <jobordner@ucsd.edu>
//
//        Date: 2007-04-10
//
//======================================================================

#include <mpi.h>
#include <assert.h>
#include <string>
#include <vector>
#include <map>

#include "HYPRE_sstruct_ls.h"

#include "mpi.hpp"
#include "scalar.hpp"
#include "point.hpp"
#include "discret.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"
#include "sphere.hpp"
#include "problem.hpp"
#include "hypre.hpp"

const int debug = 1;
const int debug_hypre = 1;
const int debug_discret = 1;

//----------------------------------------------------------------------

Hypre::Hypre (Mpi & m)
  : mpi(m)
{
  
}

//----------------------------------------------------------------------

void Hypre::init_hierarchy (Hierarchy & hierarchy)
// Initialize the grid hierarchy
{

  //======================================================================
  // CREATE GRIDS
  //======================================================================

  if (debug) printf ("DEBUG ==================\n");
  if (debug) printf ("DEBUG %s:%d CREATE GRIDS\n",__FILE__,__LINE__);
  if (debug) printf ("DEBUG ==================\n");

  for (int i=0; i<hierarchy.num_grids(); i++) {

    Grid & grid = hierarchy.grid(i);

    if (grid.is_local(mpi.ip())) {

      if (debug) printf ("DEBUG %s:%d grid %d\n",__FILE__,__LINE__,i);
      // Create the HYPRE grid structure

      if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructGridCreate (\n");
      if (debug_hypre) printf ("DEBUG_HYPRE      hierarchy.dimension() = %d\n", 
			       hierarchy.dimension());
      if (debug_hypre) printf ("DEBUG_HYPRE      hierarchy.num_grids() = %d\n", 
			       hierarchy.num_grids());

      HYPRE_SStructGridCreate (MPI_COMM_WORLD, 
			       hierarchy.dimension(), 
			       hierarchy.num_levels()+1, 
			       &grid.hypre_grid());
    }

  }

  //======================================================================
  // SET GRID EXTENTS
  //======================================================================

  if (debug) printf ("DEBUG ======================\n");
  if (debug) printf ("DEBUG %s:%d SET GRID EXTENTS\n",__FILE__,__LINE__);
  if (debug) printf ("DEBUG ======================\n");

  for (int i=0; i<hierarchy.num_grids(); i++) {

    Grid & grid = hierarchy.grid(i);

    if (grid.is_local(mpi.ip())) {

      // Set HYPRE grid extents

      int lower_extents[3] = { 0, 0, 0};
      int upper_extents[3] = {grid.num_unknowns(0)-1,
			      grid.num_unknowns(1)-1,
			      grid.num_unknowns(2)-1};

      if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructGridSetExtents\n");
      if (debug_hypre) printf ("DEBUG_HYPRE    grid.id()     = %d\n",grid.id());
      if (debug_hypre) printf ("DEBUG_HYPRE   lower_extents  = %d %d %d\n",
			       lower_extents[0],lower_extents[1],lower_extents[2]);
      if (debug_hypre) printf ("DEBUG_HYPRE    upper_extents = %d %d %d\n",
			       upper_extents[0],upper_extents[1],upper_extents[2]);

      HYPRE_SStructGridSetExtents(grid.hypre_grid(),
				  grid.level(),
				  lower_extents,
				  upper_extents);
    }

  }

  //======================================================================
  // SET GRID VARIABLE TYPE
  //======================================================================

  if (debug) printf ("DEBUG ============================\n");
  if (debug) printf ("DEBUG %s:%d SET GRID VARIABLE TYPE\n",__FILE__,__LINE__);
  if (debug) printf ("DEBUG ============================\n");

  for (int i=0; i<hierarchy.num_grids(); i++) {

    Grid & grid = hierarchy.grid(i);

    // Set HYPRE grid variables

    if (grid.is_local(mpi.ip())) {

      HYPRE_SStructVariable variable_types[] = { HYPRE_SSTRUCT_VARIABLE_CELL };

      if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructGridSetVariables(\n");
      if (debug_hypre) printf ("DEBUG_HYPRE         grid.id() = %d\n",grid.id());
      if (debug_hypre) printf ("DEBUG_HYPRE                     1 \n");
      if (debug_hypre) printf ("DEBUG_HYPRE    variable_types = %d\n",variable_types[0]);

      HYPRE_SStructGridSetVariables(grid.hypre_grid(),
				    grid.level(),
				    1,
				    variable_types);
    }
  }

  //======================================================================
  // SET GRID NEIGHBORS
  //======================================================================

  if (debug) printf ("DEBUG ========================\n");
  if (debug) printf ("DEBUG %s:%d SET GRID NEIGHBORS\n",__FILE__,__LINE__);
  if (debug) printf ("DEBUG ========================\n");

  for (int i=0; i<hierarchy.num_grids(); i++) {

    Grid & grid = hierarchy.grid(i);

    // Set HYPRE grid neighbors

    if (grid.is_local(mpi.ip())) {

      int index_map[] = {0,1,2}; // local coordinate axes same between all grids

      for (int j=0; j < grid.num_neighbors(); j++) {
	Grid & neighbor = grid.neighbor(j);
	int gl[3],gu[3];
	int nl[3],nu[3];
	grid.find_neighbor_indices (neighbor,gl,gu,nl,nu);

// 	if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructGridSetNeighborBox()\n");
// 	if (debug_hypre) printf ("DEBUG_HYPRE     grid.id()         = %d\n", grid.id());
// 	if (debug_hypre) printf ("DEBUG_HYPRE     gl[0],gl[1],gl[2] = %d %d %d\n",
// 				 gl[0],gl[1],gl[2]);
// 	if (debug_hypre) printf ("DEBUG_HYPRE     gu[0],gu[1],gu[2] = %d %d %d\n",
// 				 gu[0],gu[1],gu[2]);
// 	if (debug_hypre) printf ("DEBUG_HYPRE         neighbor.id() = %d\n",neighbor.id());
// 	if (debug_hypre) printf ("DEBUG_HYPRE     nl[0],nl[1],nl[2] = %d %d %d\n",
// 				 nl[0],nl[1],nl[2]);
// 	if (debug_hypre) printf ("DEBUG_HYPRE     nu[0],nu[1],nu[2] = %d %d %d\n",
// 				 nu[0],nu[1],nu[2]);
// 	if (debug_hypre) printf ("DEBUG_HYPRE             index_map = %d %d %d\n",index_map[0],index_map[1],index_map[2]);

	//---------------------------------------
	// Update grid's Discret's boundary mask
	//---------------------------------------

	// Find axis and face

	int axis, face;
	if (gl[0] == -1)                   {axis = 0; face = 0;}
	if (gl[0] == grid.num_unknowns(0)) {axis = 0; face = 1;}
	if (gl[1] == -1)                   {axis = 1; face = 0;}
	if (gl[1] == grid.num_unknowns(1)) {axis = 1; face = 1;}
	if (gl[2] == -1)                   {axis = 2; face = 0;}
	if (gl[2] == grid.num_unknowns(2)) {axis = 2; face = 1;}

	// Determine index bounds 

	int in0, in1, jn0, jn1;
	in0 = gl[(axis+1)%3];
	in1 = gu[(axis+1)%3];
	jn0 = gl[(axis+2)%3];
	jn1 = gu[(axis+2)%3];

	if (debug_discret) printf ("DEBUG_DISCRET axis=%d face=%d [(%d,%d) - (%d,%d)]\n",
				   axis,face,in0,jn0,in1,jn1);
      
	// Set Discret boundary mask

	for (int in=in0; in<=in1; in++) {
	  for (int jn=jn0; jn<=jn1; jn++) {
	    grid.discret().neighbor_cell(axis,face,in,jn) = Discret::_same_;
	  }
	}

	if (debug_discret) grid.discret().print();

	//---------------------------------------
	// Update HYPRE grid with neighboring box
	//---------------------------------------

	//	HYPRE_SStructGridSetNeighborBox(grid.hypre_grid(),
	//					grid.id(),    gl,gu,
	//					neighbor.id(),nl,nu,
	//					index_map);
      }
    }
  }

  //======================================================================
  // ASSEMBLE GRIDS
  //======================================================================

  if (debug) printf ("DEBUG ====================\n");
  if (debug) printf ("DEBUG %s:%d ASSEMBLE GRIDS\n",__FILE__,__LINE__);
  if (debug) printf ("DEBUG ====================\n");

  for (int i=0; i<hierarchy.num_grids(); i++) {

    Grid & grid = hierarchy.grid(i);

    if (grid.is_local(mpi.ip())) {

      // Assemble grid
      if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructGridAssemble()\n",
			       __FILE__,__LINE__);
      HYPRE_SStructGridAssemble (grid.hypre_grid());
    }
  }
  if (debug) printf ("DEBUG ==================================\n");
  if (debug) printf ("DEBUG %s:%d EXIT Hypre::init_hierarchy()\n",__FILE__,__LINE__);
  if (debug) printf ("DEBUG ==================================\n");
}

// --------------------------------------------------
// Initialize the stencils
// --------------------------------------------------

void Hypre::init_stencil (Hierarchy & hierarchy)
// Initialize the discretization stencil
{

  int dim = hierarchy.dimension();

  if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructStencilCreate()\n",
			   __FILE__,__LINE__);
  HYPRE_SStructStencilCreate (dim,dim*2+1,&stencil_);

  int x0[] = { 0,0,0 };
  if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructStencilSetEntry()\n",
			   __FILE__,__LINE__);
  HYPRE_SStructStencilSetEntry (stencil_, 0, x0, 0);

  if (dim >= 1) {
    int xp[] = {  1,0,0 };
    int xm[] = { -1,0,0 };
    if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructStencilSetEntry()\n",
			     __FILE__,__LINE__);
    HYPRE_SStructStencilSetEntry (stencil_, 1, xp, 0);
    if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructStencilSetEntry()\n",
			     __FILE__,__LINE__);
    HYPRE_SStructStencilSetEntry (stencil_, 2, xm, 0);
  }

  if (dim >= 2) {
    int yp[] = { 0, 1,0 };
    int ym[] = { 0,-1,0 };
    if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructStencilSetEntry()\n",
			     __FILE__,__LINE__);
    HYPRE_SStructStencilSetEntry (stencil_, 3, yp, 0);
    if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructStencilSetEntry()\n",
			     __FILE__,__LINE__);
    HYPRE_SStructStencilSetEntry (stencil_, 4, ym, 0);
  }

  if (dim >= 3) {
    int zp[] = { 0,0, 1 };
    int zm[] = { 0,0,-1 };
    if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructStencilSetEntry()\n",
			     __FILE__,__LINE__);
    HYPRE_SStructStencilSetEntry (stencil_, 5, zp, 0);
    if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructStencilSetEntry()\n",
			     __FILE__,__LINE__);
    HYPRE_SStructStencilSetEntry (stencil_, 6, zm, 0);
  }

}

// --------------------------------------------------
// Initialize the graphs
// --------------------------------------------------

void Hypre::init_graph (Hierarchy & hierarchy)
{
  int i;

  for (i=0; i<hierarchy.num_grids(); i++) {

    Grid & grid = hierarchy.grid(i);

    if (grid.is_local(mpi.ip())) {

      // Create the hypre graph for the grid

      if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructGraphCreate()\n",
			       __FILE__,__LINE__);
      HYPRE_SStructGraphCreate (MPI_COMM_WORLD, 
				grid.hypre_grid(), 
				&grid.hypre_graph());

      // Set the stencil for the grid

      if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructGraphSetStencil()\n",
			       __FILE__,__LINE__);
      HYPRE_SStructGraphSetStencil (grid.hypre_graph(),
				    grid.id(),
				    0,
				    stencil_);

      // Set any additional entries between grid and parent

      if (grid.level() > 0) {

      
	//      HYPRE_SStructGraphAddEntries (grid.hypre_graph(),
	//				    grid.id());

      }

      // Assemble the graph

      if (debug_hypre) printf ("DEBUG_HYPRE HYPRE_SStructGraphAssemble()\n",
			       __FILE__,__LINE__);
      HYPRE_SStructGraphAssemble (grid.hypre_graph());
    }
    
				  
  }
}

// --------------------------------------------------
// Initialize the matrix A
// --------------------------------------------------

void Hypre::init_matrix (Hierarchy & hierarchy)
{
}

// --------------------------------------------------
// Initialize the right-hand-side vector b
// --------------------------------------------------

void Hypre::init_rhs (Hierarchy & hierarchy)
{
}

// --------------------------------------------------
// Initialize the linear solver
// --------------------------------------------------

void Hypre::init_solver (Hierarchy & hierarchy)
{
}

// --------------------------------------------------
// Solve the linear system Ax = b
// --------------------------------------------------

void Hypre::solve (Hierarchy & hierarchy)
{
}

