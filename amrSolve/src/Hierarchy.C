//======================================================================
//
// File:        Hierarchy.C
//
// Package:     amrSolve
//
// Description: Interface to Enzo hierarchies
//
// Classes:     Hierarchy
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
#include "Hierarchy.h"  
#include "ItGrids.h"
#include "ItLevels.h"

  //======================================================================
  // Hierarchy
  //======================================================================

  void Hierarchy::attach (LevelHierarchyEntry * hierarchy []) 
  {

    // For each enzo-level in the enzo-hierarchy...

    for (int i=0; hierarchy[i]; i++) {

      // Create a new Level, ...

      Level * l = new Level () ;

      // ...attach it to the enzo-level, ...

      l->attach(hierarchy[i]) ;

      // ...add it to the list of Levels in the Hierarchy, ...

      levels_.push_back(l) ;

      // ...set the Grid's hierarchy to this Hierarchy

      printf ("%s:%d: LEFT OFF HERE--NEED GRID NOT LEVEL\n",
	       __FILE__,__LINE__);
      //      assert (l->phierarchy_ == 0); 
      //      l->phierarchy_ = this;

      // ...find the grid's parent

      printf ("%s:%d oops: Left off here\n",__FILE__,__LINE__);
    }
  }

  //----------------------------------------------------------------------

  void Hierarchy::detach ()
  {
    for (int i=0 ; i<levels_.size() ; i++) {
      levels_[i]->detach();
      delete levels_[i] ;
    }
  }

  //----------------------------------------------------------------------

  void Hierarchy::geomview (FILE *fpr, bool full) const 
  {
    // Color mapping for levels

    int bcolor[] = {1, 1, 0, 0, 0, 1, 1};
    int rcolor[] = {1, 0, 1, 0, 1, 0, 0};
    int gcolor[] = {1, 0, 0, 1, 1, 1, 1};

    if (full) { // ADD GEOMVIEW OBJECT?  SHARED CODE IN LEVEL.C

      // Determine number of grids in the Hierarchy

      ItLevels itL (*this);
      int n = 0;
      while (const Level * l = ++itL) {
	n += l->numGrids();
      }

      // Print first two lines of geomview *.vect file
      fprintf (fpr,"VECT\n");
      fprintf (fpr,"%d %d %d\n",4*n, 16*n, this->numLevels());

      //
      for (int i=0; i<n; i++) fprintf (fpr,"8 3 3 2 "); fprintf (fpr,"\n");

      itL.reset();
      while (const Level *l = ++itL) {
	fprintf (fpr,"1 0 0 0 ");
	for (int i=1; i<l->numGrids(); i++) {
	  fprintf (fpr,"0 0 0 0 "); 
	}
	fprintf (fpr,"\n");
      }
    }

    ItLevels itL (*this);
    while (const Level * l = ++itL) {
      ItGrids itG (*l);
      while (const Grid * g = ++itG) {
	g->geomview(fpr,0);
      }
    }

    if (full) {
      for (int i=0; i<this->numLevels(); i++) {
	int j=i%7;
	fprintf (fpr,"%d %d %d 0\n",rcolor[j],gcolor[j],bcolor[j]);
      }
    }
  }

} // namespace
