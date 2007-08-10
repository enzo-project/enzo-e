
/// Hierarchy class source file

/**
 * 
 * 
 * 
 * @file
 * @author James Bordner <jobordner@ucsd.edu>
 *
 */

#include <assert.h>
#include <stdio.h>
#include <vector>
#include <map>

#include "HYPRE_sstruct_ls.h"

#include "scalar.hpp"
#include "discret.hpp"
#include "mpi.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"

//----------------------------------------------------------------------

const int debug = 1;

//----------------------------------------------------------------------

Hierarchy::Hierarchy () throw ()
  : dimension_(0)
{
}
	  
//----------------------------------------------------------------------

Hierarchy::~Hierarchy () throw ()
{
}

//======================================================================

void Hierarchy::insert_grid (Grid * pgrid) throw ()
{
  printf ("Hierarchy::insert_grid()\n");
  if (pgrid->id() >= grids_.size()) {
    if (debug) printf ("DEBUG: resizing Hierarchy::grids_ from %d to %d\n",
		       grids_.size(),pgrid->id() + 1);
    grids_.resize (pgrid->id() + 1);
  }
  grids_[pgrid->id()] = pgrid;

}

//----------------------------------------------------------------------

void Hierarchy::init_levels () throw ()
{
  printf ("Hierarchy::init_levels()\n");

  //-------------------------
  //  Initialize grid parents
  //-------------------------

  int i,j,j1,j2,k;

  for (i=0; i<num_grids (); i++) {
    Grid * g = &grid(i);
    Grid * p = (g->id_parent() >= 0) ? & grid(g->id_parent()) : 0;
    this->set_parent(g,p);
  }

  //------------------------
  //  Initialize grid levels
  //------------------------

  bool done = false;
  while (! done) {
    done = true;
    for (i=0; i<num_grids (); i++) {
      Grid * g = &grid(i);
      if (g->level() < 0) { // Haven't determined level yet
	if (parent(g) == 0) {
	  // Grids without parents are in level 0
	  g->set_level(0);
	  insert_in_level_ (0,*g);
	} else if (parent(g)->level() >= 0) {
	  // Grids with parents of known level have level = 1 + parent level
	  int level = parent(g)->level() + 1;
	  g->set_level(level);
	  insert_in_level_ (level,*g);
	} else {
	  // Otherwise grids have parents of unknown level
	  done = false;
	}
      }
      g->print();
    }
  }

  //--------------------------
  //  Initialize grid children
  //--------------------------

  for (i=0; i<num_grids(); i++) {
    Grid * g = & grid(i);
    // If a grid has a parent, then the grid is the parent's child
    if (parent(g) != 0) {
      parent(g)->set_child(*g);
    }
  }

  //---------------------------
  //  Initialize grid neighbors
  //---------------------------

  // First level 0: test all pairs

  for (i=0; i<level(0).num_grids(); i++) {
    Grid * g1 = & level(0).grid(i);
    for (j=i+1; j<level(0).num_grids(); j++) {
      Grid * g2 = & level(0).grid(j);
      if (g1->is_adjacent(*g2)) {
	if (debug) printf ("DEBUG grids %d and %d are adjacent\n",
			   g1->id(),g2->id());
	g1->set_neighbor (*g2);
	g2->set_neighbor (*g1);
      }
    }
  }

  // Next levels > 0: for each grid, only test against
  // parents' children, and parents' neighbors' children,
  // since if two grids' parents are not adjacent, then
  // they necessarily cannot be neighbors

  for (k=1; k<num_levels(); k++) {

    for (i=0; i<level(k).num_grids(); i++) {

      Grid * g1 = & level(k).grid(i);

      // Check parents' children

      for (j=0; j<parent(g1)->num_children(); j++) {
	Grid * g2 = & parent(g1)->child(j);
	if (g1->is_adjacent(*g2) && g1->id() > g2->id()) {
	  assert_neighbors (*g1,*g2);
	}
      }

      // Check parents' neighbors' children

      for (j1=0; j1<parent(g1)->num_neighbors(); j1++) {
	Grid * gn = & parent(g1)->neighbor(j1);
	for (j2=0; j2<gn->num_children(); j2++) {
	  Grid * g2 = & gn->child(j2);
	  if (g1->is_adjacent(*g2) && g1->id() > g2->id()) {
	    if (debug) printf ("DEBUG grids %d and %d are adjacent cousins\n",
			       g1->id(),g2->id());
	    assert_neighbors (*g1,*g2);
	  }
	}
      }
    }
  }
  //   printf ("%s:%d Grid::init_mesh()\n",__FILE__,__LINE__);

}

//======================================================================

void Hierarchy::print () throw ()
{
  printf ("Hierarchy\n");
  for (int i=0; i<num_levels(); i++) {
    level(i).print();
  }
}

//----------------------------------------------------------------------

void Hierarchy::write (FILE *fp) throw ()
{
  if (fp == 0) fp = stdout;

  fprintf (fp,"Hierarchy\n");
  for (int i=0; i<num_levels(); i++) {
    level(i).write();
  }
}

//======================================================================
// PRIVATE MEMBER FUNCTIONS
//======================================================================

void Hierarchy::insert_in_level_ (int level, Grid & grid) throw ()
{
  // Resize levels_[] if needed
  if (level >= levels_.size()) {
    if (debug) printf ("DEBUG: resizing Hierarchy::levels_ from %d to %d\n",
		       levels_.size(),level + 1);
    levels_.resize (level + 1);
  }
  if (levels_.at(level) == 0) {
    if (debug) printf ("DEBUG: creating new Level at %d\n",level);
    levels_[level] = new Level(level);
  }
  levels_[level]->insert_grid (grid);
}
