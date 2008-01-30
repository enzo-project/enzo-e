//345678901234567890123456789012345678901234567890123456789012345678901234567890

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
#include <math.h>

#include <vector>
#include <map>

#include "HYPRE_sstruct_ls.h"

#include "hypre-solve.hpp"

#include "scalar.hpp"
#include "error.hpp"
#include "point.hpp"
#include "domain.hpp"
#include "faces.hpp"
#include "mpi.hpp"
#include "grid.hpp"
#include "level.hpp"
#include "hierarchy.hpp"


//----------------------------------------------------------------------

const int trace          = 0;
const int debug          = 0;
const int debug_detailed = 0;
const int geomview       = 1;

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
			    Mpi    & mpi) throw ()
{
  if (debug) printf ("Hierarchy::init_levels()\n");

  init_grid_parents_();
  init_grid_levels_();
  init_grid_children_();
  init_grid_neighbors_();
  init_indices_();                // DEPENDENCY: Requires init_grid_levels_()
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

      if (g1->is_adjacent(*g2)) {
	g1->set_neighbor (*g2);
	g2->set_neighbor (*g1);
	if (debug) printf ("DEBUG grids %d and %d are adjacent\n",
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
	if (g1->is_adjacent(*g2) && g1->id() > g2->id()) {
	  assert_neighbors (*g1,*g2);
	}
      }

      // Check parents' neighbors' children

      for (j1=0; j1<parent(*g1)->num_neighbors(); j1++) {
	Grid * gn = & parent(*g1)->neighbor(j1);
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

  while (Level * level = itl++) {
    ItLevelGridsAll itg (*level);

    while (Grid * grid = itg++) {
      ItGridNeighbors itn (*grid);

      // Set "adjacent" pointers for adjacent grids in same level

      Grid * adjacent;

      while (adjacent = itn++) {
	if (grid->is_local() || adjacent->is_local()) {
	  if (grid->neighbor_shared_face
	      (*adjacent,axis,face,il0,il1,iu0,iu1)) {
	    for (ig0=il0; ig0<iu0; ig0++) {
	      for (ig1=il1; ig1<iu1; ig1++) {
		grid->faces().adjacent(axis,face,ig0,ig1) = adjacent;
	      }
	    }
	  }
	}
      }

      // Set "adjacent" pointers for adjacent non-parent coarse grids

      if (parent (*grid)) {

	ItGridNeighbors itpn (*parent(*grid));

	while (adjacent = itpn++) {
	  _TRACE_;
	  if (grid->is_local() || adjacent->is_local()) {
	    if (grid->coarse_shared_face
		(*adjacent,axis,face,il0,il1,iu0,iu1)) {
	      _TRACE_;
	      if (debug) printf ("%s:%d %d %d  %d %d\n",__FILE__,__LINE__,il0,iu0,il1,iu1);
	      for (ig0=il0; ig0<iu0; ig0++) {
		for (ig1=il1; ig1<iu1; ig1++) {
		  grid->faces().adjacent(axis,face,ig0,ig1) = adjacent;
		}
	      }
	    }
	  }
	}

	// Set remaining "adjacent" pointers to parent grid

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
		if (grid->faces().adjacent(axis,face,ig0,ig1) == NULL) {
		  grid->faces().adjacent(axis,face,ig0,ig1) = parent (*grid);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // ------------------------------------------------------------
  // 1. Categorize boundary face-zones
  // ------------------------------------------------------------

  int ih0[3][2];
  this->indices0(ih0);
  ItLevelGridsLocal itgl (this->level(0));
  while (Grid * grid = itgl++) {
    int ig[3][2];
    grid->indices(ig);
    for (axis = 0; axis < 3; axis++) {
      for (face = 0; face < 2; face++) {
	if (ih0[axis][face] == ig[axis][face]) {
	  grid->faces().label(axis,face,Faces::_boundary_);
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
    ItLevelGridsAll itg (*level);
    while (Grid * grid = itg++) {
      //      int ig3[3][2];
      //      grid->indices(ig3);
      ItGridNeighbors itn (*grid);
      while (Grid * neighbor = itn++) {
	if (grid->neighbor_shared_face 
	    (*neighbor, axis, face, il0,il1,iu0,iu1)) {
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

      //      if (mpi.is_root()) {
      //	sprintf (filename,"grid-L%d.vect",level->index());
      //	FILE * fp = fopen (filename,"w");
      //	sprintf (filename,"grid-L%d.vect",level->index());
      //	level->geomview_grid (fp);
      //      }
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

/// Write the Hierarchy grids to the given open file in geomview format

// void Hierarchy::geomview_grid_ (FILE *fpr, bool full) throw ()
// {

//   int ng = num_grids();
//   int nl = num_levels();

//   // Write header

//   if (full) {

//     // Print first two lines of geomview *.vect file
//     fprintf (fpr,"VECT\n");
//     fprintf (fpr,"%d %d %d\n",4*ng, 16*ng, nl);

//     //
//     for (int i=0; i<ng; i++) fprintf (fpr,"8 3 3 2 "); fprintf (fpr,"\n");

//     ItHierarchyLevels itl(*this);
//     while (Level *level = itl++) {
//       fprintf (fpr,"1 0 0 0 ");
//       for (int i=1; i<level->num_grids(); i++) {
// 	fprintf (fpr,"0 0 0 0 "); 
//       }
//       fprintf (fpr,"\n");
//     }
//   }

//   // For each level, print out all grids in the level

//   ItHierarchyLevels itl(*this);

//   while (Level * level = itl++) {

//     ItLevelGridsAll itg (*level);

//     while (Grid * grid = itg++) {

//       grid->geomview_grid(fpr,0);

//     }

//   }

//   // Color mapping for levels

//   int bcolor[] = {1, 1, 0, 0, 0, 1, 1};
//   int rcolor[] = {1, 0, 1, 0, 1, 0, 1};
//   int gcolor[] = {1, 0, 0, 1, 1, 1, 0};

//   if (full) {
//     for (int i=0; i<nl; i++) {
//       int j=i%7;  // 7 is length of [rgb]color[] arrays
//       fprintf (fpr,"%d %d %d 0\n",rcolor[j],gcolor[j],bcolor[j]);
//     }
//   }
// }

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
