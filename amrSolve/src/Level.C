//======================================================================
//
// File:        Level.C
//
// Package:     amrSolve
//
// Description: Interface to Enzo levels
//
// Classes:     Level
//
//----------------------------------------------------------------------
//
// Author:      James Bordner (jbordner@cosmos.ucsd.edu)
//
// History:     2002-03-31  Created
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//======================================================================

#include <assert.h>

#include "enzo-grid.h"
#include "enzo-level.h"

namespace amrSolve {

#include "Level.h"
#include "ItGrids.h"

  //======================================================================
  // Level
  //======================================================================

  void Level::attach (LevelHierarchyEntry * enzo_level) 
  {
    // Add Grids to Level

    for (; enzo_level; enzo_level = enzo_level -> NextGridThisLevel) {

      // Create a new Grid, ...

      Grid * g = new Grid ();

      // ...attach it to the enzo-grid, ...

      g->attach(enzo_level->GridData);

      // ...add it to the list of Grids in the Level, ...

      grids_.push_back(g);

      // ...set the Grid's level to this Level

      assert (g->plevel_ == 0); 
      g->plevel_ = this;
      
    }

    // Find neighboring Grids in level

    ItGrids itG1 (*this);
    while (const Grid * g1 = ++itG1) {
      ItGrids itG2 (*this);
      while (const Grid * g2 = ++itG2) {
	if (g1->isNeighbor(*g2)) {
	  assertNeighbors (*g1,*g2);
	}
      }
    }
  }

  //----------------------------------------------------------------------

  void Level::detach ()
  {
    for (int i=0; i<grids_.size(); i++) {
      grids_[i]->detach();
      delete grids_[i];
    }
  }

  //----------------------------------------------------------------------

  void Level::geomview (FILE *fpr, bool full) const 
  {
    int n = this->numGrids();

    if (full) { // ADD GEOMVIEW OBJECT?  SHARED CODE IN LEVEL.C
      fprintf (fpr,"VECT\n");
      fprintf (fpr,"%d %d 1\n",4*n, 16*n);
      for (int i=0; i<n; i++) fprintf (fpr,"8 3 3 2 "); fprintf (fpr,"\n");
      fprintf (fpr,"1 0 0 0 ");
      for (int i=1; i<n; i++) fprintf (fpr,"0 0 0 0 "); fprintf (fpr,"\n");
    }

    ItGrids itG (*this);
    while (const Grid *g = ++itG) {
      g->geomview(fpr,0);
    }

    if (full) {
      fprintf (fpr,"1 1 1 0\n");
    }
  }

} // namespace
