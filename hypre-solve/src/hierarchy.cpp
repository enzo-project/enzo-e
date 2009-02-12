//345678901234567890123456789012345678901234567890123456789012345678901234567890

/// Hierarchy class source file

/**
 * 
 * @file      hierarchy.cpp
 * @brief     Implementation of the Hierarchy class
 * @author    James Bordner
 * @bug       none
 *
 * $Id$
 *
 */

#include <assert.h>
#include <stdio.h>
#include <math.h>

#include <vector>
#include <map>
#include <string>

#include "HYPRE_sstruct_ls.h"

#include "newgrav-hypre-solve.h"

//----------------------------------------------------------------------

const int trace          = 0;
const int debug          = 0;
const int debug_detailed = 0;
const int geomview       = 0;

//----------------------------------------------------------------------

#include "newgrav-scalar.h"
#include "newgrav-error.h"
#include "newgrav-point.h"
#include "newgrav-domain.h"
#include "newgrav-faces.h"
#include "newgrav-mpi.h"
#include "newgrav-grid.h"
#include "newgrav-level.h"
#include "newgrav-hierarchy.h"


//----------------------------------------------------------------------

/// Hierarchy class constructor.

/** Currently does nothing */

Hierarchy::Hierarchy () throw ()
  : dimension_(0)
{
  grids0_.push_back(0);
  levels0_.push_back(0);
}
	  
//----------------------------------------------------------------------

Hierarchy::~Hierarchy () throw ()

/// Hierarchy class descructor.

/** Currently does nothing */

{
}

//======================================================================

void Hierarchy::enzo_attach () throw ()
{
  TEMPORARY("Hierarchy::enzo_attach()");
  set_dim(3);
  // For each grid...
  //
  //  create "value" string consisting of the following:
  //
  //	 id:                Unique id, starting from 0 and counting up
  //     id_parent:         id of the parent
  //     ip                 MPI processor rank
  //     xl[0] xl[1] xl[2]  Lower position defining the grid extent
  //     xu[0] xu[1] xu[2]  Upper position defining the grid extent
  //     il[0] il[1] il[2]  Coordinate defining the lower grid index
  //     n[0] n[1] n[2]     Size of the grid
  //
  //  insert the grid into the hierarchy:
  //
  //     hierarchy_.insert_grid(new Grid(value));
  
}

//----------------------------------------------------------------------

void Hierarchy::enzo_detach () throw ()
{
  TEMPORARY("Hierarchy::enzo_detach()");
}

//======================================================================

/// Insert a grid into the hierarchy.

/** Inserts the grid's pointer into the grids0_ vector, at the
  * position given by the grid's id(), and always maintaining a 0
  * sentinel at the end of grids0_ */

void Hierarchy::insert_grid (Grid * pgrid) throw ()
{
  // Insert grid into the list of all grids in the Hierarchy

  if (unsigned(pgrid->id()) + 1 >= grids0_.size()) {
    grids0_.resize (pgrid->id() + 2);
    grids0_[grids0_.size() - 1] = 0;
  }
  grids0_[pgrid->id()] = pgrid;
}

//======================================================================

/// Initialize the hierarchy after all grids have been inserted

/** After all grids are inserted into the hierarchy, this function
    determines each grid's parent, containing level, children, 
    neighbors, face-zone categories, and global indices. */

void Hierarchy::initialize (Domain & domain,
			    Mpi    & mpi,
			    bool   is_periodic) throw ()
{
  if (debug) printf ("Hierarchy::init_levels()\n");

  init_grid_parents_();
  init_grid_levels_();
  init_grid_children_();
  init_grid_neighbors_();
  init_indices_();                // DEPENDENCY: Requires init_grid_levels_()
  init_extents_(is_periodic);     // DEPENDENCY: Requires init_grid_levels_()
  init_grid_faces_(domain, mpi);  // DEPENDENCY: Requires init_indices_()

  geomview_grids(mpi);

}

//------------------------------------------------------------------------

/// Determines the parent grid of each grid in the hierarchy.

void Hierarchy::init_grid_parents_ () throw ()
{
  ItHierarchyGridsAll itg (*this);
  while (Grid * g = itg++) {
    Grid * p = (g->id_parent() >= 0) ? & grid(g->id_parent()) : 0;
    this->set_parent(g,p);
  }
}

//------------------------------------------------------------------------

/// Determines the level of each grid in the hierarchy.

void Hierarchy::init_grid_levels_ () throw ()
{
  bool done = false;

  while (! done) {

    done = true;

    // Loop through grids, and try to determine their level if unknown

    ItHierarchyGridsAll itg (*this);
    while (Grid * g = itg++) {

      // If grid's level < 0, then we haven't determined its level yet

      if (g->level() < 0) { 

	// Grids without parents must be in level 0 ...

	if (parent(*g) == 0) {
	  g->set_level(0);
	  insert_in_level_ (0,*g);
	}
	
	// ... grids with parents of known level have level = 1 +
	// parent level ...

	else if (parent(*g)->level() >= 0) {
	  int level = parent(*g)->level() + 1;
	  g->set_level(level);
	  insert_in_level_ (level,*g);
	} 
	
	// ... otherwise a grid's parents is in an unknown level, so
	// we're not done yet

	else {
	  done = false;
	}
      }
    }
  }
}

//------------------------------------------------------------------------

///  Initialize grid children.

void Hierarchy::init_grid_children_ () throw ()
{
  ItHierarchyGridsAll itg (*this);
  while (Grid * g = itg++) {
    // If a grid has a parent, then the grid is the parent's child
    if (parent(*g) != 0) {
      parent(*g)->set_child(*g);
    }
  }
}


//------------------------------------------------------------------------

///  Find each grid's neighbors.

/**  If a grid is in level 0, then find its neighbors by comparing
     with all other grids in the level.  Otherwise, if a grid is in a
     level > 0, then we can save time by only testing grids that are
     children of the parent, and children of all the parents
     neighbors. */

void Hierarchy::init_grid_neighbors_ () throw ()
{

  int i,j;

  // For level == 0, test all pairs

  for (i=0; i<level(0).num_grids(); i++) {

    Grid * g1 = & level(0).grid(i);

    for (j=i+1; j<level(0).num_grids(); j++) {

      Grid * g2 = & level(0).grid(j);

      if (g1->is_adjacent(*g2,period_)) {
	assert_neighbors(*g1,*g2);
	if (debug) printf ("DEBUG 1 grids %d and %d are neighbors\n",
			   g1->id(),g2->id());
      }
    }
  }

  // For levels > 0, only test parents' children, and parents'
  // neighbors' children.  If two grids' parents are not the same or
  // not neighbors, then they necessarily cannot be neighbors.

  int k,j1,j2;

  for (k=1; k<num_levels(); k++) {

    for (i=0; i<level(k).num_grids(); i++) {

      Grid * g1 = & level(k).grid(i);

      // Check parents' children

      for (j=0; j<parent(*g1)->num_children(); j++) {
	Grid * g2 = & parent(*g1)->child(j);
	if (g1->is_adjacent(*g2,period_) && g1->id() > g2->id()) {
	  if (debug) printf ("DEBUG 2 grids %d and %d are neighbors\n",
			     g1->id(),g2->id());
	  assert_neighbors (*g1,*g2);
	}
      }

      // Check parents' neighbors' children

      for (j1=0; j1<parent(*g1)->num_neighbors(); j1++) {
	Grid * gn = & parent(*g1)->neighbor(j1);
	for (j2=0; j2<gn->num_children(); j2++) {
	  Grid * g2 = & gn->child(j2);
	  if (g1->is_adjacent(*g2,period_) && g1->id() > g2->id()) {
	    if (debug) printf ("DEBUG 3 grids %d and %d are neighbors\n",
			       g1->id(),g2->id());
	    assert_neighbors (*g1,*g2);
	  }
	}
      }
    }
  }
}

//------------------------------------------------------------------------

///  Find each grid's face categories.

/** After all grid inter-connections are determined, this function
    determines the neighbor structure for each face-zone. Only
    performed for local grids, and grids neighboring local grids. 
    Allocates Faces for grids that haven't had them allocated yet. */

void Hierarchy::init_grid_faces_ (Domain & domain,
				  Mpi    & mpi) throw ()

{
  int axis, face;
  int ig0,in0,il0,iu0;
  int ig1,in1,il1,iu1;
  int j0,j1;            // face axes

  // ------------------------------------------------------------
  // Determine pointers to neighboring grid for each face-zone
  // ------------------------------------------------------------

  ItHierarchyLevels itl (*this);
  Level * level;

  while ((level = itl++)) {

    ItLevelGridsAll itg (*level);

    // Determine level periodicity

    int period[3];
    for (int i=0; i<3; i++) {
      period[i] = period_[level->index()]==0 ? 0 : level->zones(i);
    }

    if (debug) printf ("Level %d\n",level->index());
    while (Grid * grid = itg++) {
      if (debug) grid->print();
      ItGridNeighbors itn (*grid);

      // Set "adjacent" pointers for adjacent grids in same level

      Grid * adjacent;

      while ((adjacent = itn++)) {
	if (grid->is_local() || adjacent->is_local()) {
	  if (grid->neighbor_shared_face
	      (*adjacent,axis,face,il0,il1,iu0,iu1,period)) {
	    for (ig0=il0; ig0<=iu0; ig0++) {
	      for (ig1=il1; ig1<=iu1; ig1++) {
		grid->faces().adjacent(axis,face,ig0,ig1) = adjacent;
	      }
	    }
	  }
	}
      }

      // Set "adjacent" pointers for adjacent non-parent coarse grids

      if (parent (*grid) != NULL) {

	ItGridNeighbors itpn (*parent(*grid));

	while ((adjacent = itpn++)) {
	  if (grid->is_local() || adjacent->is_local()) {
	    if (grid->coarse_shared_face
		(*adjacent,axis,face,il0,il1,iu0,iu1,period)) {
	      for (ig0=il0; ig0<=iu0; ig0++) {
		for (ig1=il1; ig1<=iu1; ig1++) {
		  grid->faces().adjacent(axis,face,ig0,ig1) = adjacent;
		}
	      }
	    }
	  }
	}

	// Set remaining "adjacent" pointers to parent grid

	if (grid->is_local() || parent(*grid)->is_local()) {
	  int num=0;

	  while (grid->parent_interior_face
		 (*parent(*grid),axis,face,il0,il1,iu0,iu1,num)) {
	    for (ig0=il0; ig0<=iu0; ig0++) {
	      for (ig1=il1; ig1<=iu1; ig1++) {
		if (grid->faces().adjacent(axis,face,ig0,ig1) == NULL) {
		  grid->faces().adjacent(axis,face,ig0,ig1) = parent (*grid);
		}
	      }
	    }
	  }
	}
      }

      if (debug_detailed) {
	int ig3[3][2];
	grid->indices(ig3);
	for (axis=0; axis<3; axis++) {
	  int j0 = (axis+1)%3;
	  int j1 = (axis+2)%3;
	  int n0 = ig3[j0][1] - ig3[j0][0];
	  int n1 = ig3[j1][1] - ig3[j1][0];
	  for (face=0; face<2; face++) {
	    for (ig0=0; ig0<n0; ig0++) {
	      for (ig1=0; ig1<n1; ig1++) {
		Grid * neighbor = grid->faces().adjacent(axis,face,ig0,ig1);
		printf ("ip=%d grid=%d neighbor=%d (axis=%d face=%d) [%d,%d]\n",
			pmpi->ip(),grid->id(),neighbor?neighbor->id():-1,
			axis,face,ig0,ig1);
	      }
	    }
	  }
	}
	fflush(stdout);
      }
      
      
    }
  }

  // ------------------------------------------------------------
  // 1. Categorize boundary face-zones
  // ------------------------------------------------------------

  // Domain boundary
  double db3[3][2]; 
  domain.lower(db3[0][0],db3[1][0],db3[2][0]);
  domain.upper(db3[0][1],db3[1][1],db3[2][1]);

  int ih0[3][2];
  this->indices0(ih0);
  while ((level = itl++)) {
    ItLevelGridsLocal itgl (*level);
    while (Grid * grid = itgl++) {

      // Grid boundary
      double gb3[3][2];
      grid->x_lower(gb3[0][0],gb3[1][0],gb3[2][0]);
      grid->x_upper(gb3[0][1],gb3[1][1],gb3[2][1]);

      // Grid spacing (should be level-independent)
      double h3[3];
      grid->h(h3[0],h3[1],h3[2]);

      for (axis = 0; axis < 3; axis++) {
	for (face = 0; face < 2; face++) {
	  if ( fabs(gb3[axis][face] - db3[axis][face]) < 0.5*h3[axis]) {
	    grid->faces().label(axis,face,Faces::_boundary_);
	  }
	}
      }
    }
  }

  // ------------------------------------------------------------
  // 2. Categorize covered (and temporary adjacent-covered) face-zones
  // ------------------------------------------------------------

  for (int ilevel = num_levels()-1; ilevel>0; ilevel--) {

    Level *level = &this->level(ilevel);

    ItLevelGridsAll itg (*level);
    while (Grid * grid = itg++) {
      Grid * parent = this->parent(*grid);
      if (grid->is_local() || parent->is_local()) {
	int num = 0;
	while (grid->parent_shared_face 
	       (*parent, axis, face, il0,il1,iu0,iu1,num)) {
	  int ig3[3][2];
	  grid->indices(ig3);
	  for (ig0=il0; ig0<=iu0; ig0++) {
	    for (ig1=il1; ig1<=iu1; ig1++) {
	      // Label covered
	      parent->faces().label(axis,face,ig0,ig1) = Faces::_covered_;
	      // Label adjacent-covered
	      Grid * neighbor = parent->faces().adjacent(axis,face,ig0,ig1);
	      if (neighbor != NULL && neighbor->level() == grid->level()) {
		int in3[3][2];
		neighbor->indices(in3);
		j0 = (axis+1)%3;
		j1 = (axis+2)%3;
		in0 = ig0 + in3[j0][0] - ig3[j0][0];
		in1 = ig1 + in3[j1][0] - ig3[j1][0];
		if (neighbor->faces().label(axis,1-face,in0,in1) != Faces::_covered_) {
		  neighbor->faces().label(axis,1-face,in0,in1) 
		    = Faces::_adjacent_covered_;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // ------------------------------------------------------------
  // 3. Categorize fine and neighbor face-zones
  // ------------------------------------------------------------

  for (int ilevel = 0; ilevel<num_levels(); ilevel++) {
    Level *level = &this->level(ilevel);

    // Determine level periodicity

    int period[3];
    for (int i=0; i<3; i++) {
      period[i] = period_[level->index()]==0 ? 0 : level->zones(i);
    }

    ItLevelGridsAll itg (*level);
    while (Grid * grid = itg++) {
      //      int ig3[3][2];
      //      grid->indices(ig3);
      ItGridNeighbors itn (*grid);
      while (Grid * neighbor = itn++) {
	if (grid->neighbor_shared_face 
	    (*neighbor, axis, face, il0,il1,iu0,iu1,period)) {
	  for (ig0=il0; ig0<=iu0; ig0++) {
	    for (ig1=il1; ig1<=iu1; ig1++) {
	      Faces::Label & fz = grid->faces().label(axis,face,ig0,ig1);
	      if (fz == Faces::_adjacent_covered_) {
		fz = Faces::_fine_;
	      } else if (fz != Faces::_covered_) {
		fz = Faces::_neighbor_;
	      }
	    }
	  }
	}
      }	  
    }
  }

  // ------------------------------------------------------------
  // 4. Categorize coarse face-zones
  // ------------------------------------------------------------

  while (Level * level = itl++) {
    ItLevelGridsAll itga(*level);
    while (Grid * grid = itga++) {
      int ig3[3][2];
      grid->indices(ig3);
      for (axis=0; axis<3; axis++) {
	int j0 = (axis+1)%3;
	int j1 = (axis+2)%3;
	int n0 = ig3[j0][1]-ig3[j0][0];
	int n1 = ig3[j1][1]-ig3[j1][0];
	for (face=0; face<2; face++) {
	  for (ig0=0; ig0<n0; ig0++) {
	    for (ig1=0; ig1<n1; ig1++) {
	      Faces::Label & fz = grid->faces().label(axis,face,ig0,ig1);
	      if (fz==Faces::_unknown_) fz = Faces::_coarse_;
	    }
	  }
	}
      }
    }
  }

  // ------------------------------------------------------------
  // Dump out any requested geomview files
  // ------------------------------------------------------------

  if (geomview) {
    char filename[80];
    for (int ilevel=0; ilevel<this->num_levels(); ilevel++) {
      for (Faces::Label type=Faces::_first_; 
	   type<=Faces::_last_; 
	   type = Faces::Label(type+1)) {
	sprintf (filename,"facezones-L%d-P%d-%s.quad",
		 ilevel,mpi.ip(),Faces::LabelName[type]);
	FILE * fp = fopen (filename,"w");
	this->level(ilevel).geomview_face_types(fp,&type,1);
	fclose(fp);
      }
    }
  }
}

//======================================================================

void Hierarchy::init_indices_ () throw()
{
  _TRACE_;
  // Determine problem size the hard way

  ItLevelGridsAll itg (level(0));

  int i,lower[3],upper[3];

  for (i=0; i<3; i++) {
    lower[i] = INT_MAX;
    upper[i] = INT_MIN;
  }

  while (Grid *grid = itg++) {
    for (i=0; i<3; i++) {
      lower[i] = MIN(grid->i_lower(i),lower[i]);
      upper[i] = MAX(grid->i_upper(i),upper[i]);
    }
  }

  for (i=0; i<3; i++) {
    il0_[i] = lower[i];
    n0_[i] = upper[i] - lower[i] + 1;
  }
  _TRACE_;
}

//======================================================================

void Hierarchy::init_extents_ (bool is_periodic) throw()
{
  _TRACE_;

  // Determine hierarchy extents from grid extents

  ItLevelGridsAll itg (level(0));

  int i;

  for (i=0; i<3; i++) {
    xl_[i] = SCALAR_MAX;
    xu_[i] = -SCALAR_MAX;
  }

  while (Grid *grid = itg++) {

    // Get grid extents

    Scalar gl[3],gu[3];
    grid->x_lower(gl[0],gl[1],gl[2]);
    grid->x_upper(gu[0],gu[1],gu[2]);

    // Adjust hierarchy extents

    for (i=0; i<3; i++) {
      xl_[i] = MIN(xl_[i],gl[i]);
      xu_[i] = MAX(xu_[i],gu[i]);
    }
  }

  for (i=0; i<3; i++) {
    period_[i] = is_periodic ? xu_[i] - xl_[i] : 0.0;
  }

  printf ("%s:%d hierarchy extents: (%g %g %g) (%g %g %g)  period: (%g %g %g)\n",
	  __FILE__,__LINE__,xl_[0],xl_[1],xl_[2],xu_[0],xu_[1],xu_[2],
	  period_[0],period_[1],period_[2]);

  _TRACE_;
}

//======================================================================

/// Write hierarchy grids to geomview files grid-L<level>-P<processor>.vect 

void Hierarchy::geomview_grids (Mpi & mpi) throw ()
{

  if (geomview) {
    ItHierarchyLevels itl (*this);
    while (Level * level = itl++) {
      char filename[20];

      sprintf (filename,"grid-L%d-P%d.vect",level->index(),mpi.ip());
      FILE * fp = fopen (filename,"w");
      level->geomview_grid_local (fp);
      fclose (fp);

    }
  }
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
  // Resize levels0_[] if needed
  if (unsigned(level + 1) >= levels0_.size()) {
    levels0_.resize (level + 2);
    levels0_[levels0_.size() - 1] = 0;
  }
  if (levels0_.at(level) == 0) {
    if (debug) printf ("DEBUG: creating new Level at %d\n",level);
    levels0_[level] = new Level(level);
  }
  levels0_[level]->insert_grid (grid);
}
