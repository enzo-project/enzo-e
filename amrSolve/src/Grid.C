//======================================================================
//
// File:        Grid.C
//
// Package:     amrSolve
//
// Description: Interface to Enzo grids
//
// Classes:     Grid
//
//----------------------------------------------------------------------
//
// Author:      James Bordner (jbordner@cosmos.ucsd.edu)
//
// History:     2002-03-22  Created
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//======================================================================


#include "enzo-grid.h"

namespace amrSolve {

#include "Grid.h"

  //======================================================================
  // Grid
  //======================================================================

  void Grid::attach (class grid * grid) 

  { grid_ = grid;

    n_[0] = grid->GridEndIndex[0] - grid->GridStartIndex[0] + 1;
    n_[1] = grid->GridEndIndex[1] - grid->GridStartIndex[1] + 1;
    n_[2] = grid->GridEndIndex[2] - grid->GridStartIndex[2] + 1;

    lp_.set(grid->GridLeftEdge[0],
	    grid->GridLeftEdge[1],
	    grid->GridLeftEdge[2]);
    hp_.set(grid->GridRightEdge[0],
	    grid->GridRightEdge[1],
	    grid->GridRightEdge[2]);
    sort (lp_, hp_); 
  }

  //----------------------------------------------------------------------

  void Grid::detach ()
  { 
    grid_ = 0;
  }

  //----------------------------------------------------------------------


  void Grid::geomview (FILE *fpr, bool full) const
  {

    if (full) {
      fprintf (fpr,"VECT\n");
      fprintf (fpr,"4 16 1\n");
      fprintf (fpr,"8 3 3 2\n");
      fprintf (fpr,"1 0 0 0\n");
    }
    fprintf (fpr,"%g %g %g\n",lp_(0),lp_(1),lp_(2));
    fprintf (fpr,"%g %g %g\n",hp_(0),lp_(1),lp_(2));
    fprintf (fpr,"%g %g %g\n",hp_(0),hp_(1),lp_(2));
    fprintf (fpr,"%g %g %g\n",lp_(0),hp_(1),lp_(2));
    fprintf (fpr,"%g %g %g\n",lp_(0),lp_(1),lp_(2));
    fprintf (fpr,"%g %g %g\n",lp_(0),lp_(1),hp_(2));
    fprintf (fpr,"%g %g %g\n",hp_(0),lp_(1),hp_(2));
    fprintf (fpr,"%g %g %g\n",hp_(0),lp_(1),lp_(2));
    fprintf (fpr,"%g %g %g\n",lp_(0),hp_(1),lp_(2));
    fprintf (fpr,"%g %g %g\n",lp_(0),hp_(1),hp_(2));
    fprintf (fpr,"%g %g %g\n",lp_(0),lp_(1),hp_(2));
    fprintf (fpr,"%g %g %g\n",hp_(0),hp_(1),lp_(2));
    fprintf (fpr,"%g %g %g\n",hp_(0),hp_(1),hp_(2));
    fprintf (fpr,"%g %g %g\n",hp_(0),lp_(1),hp_(2));
    fprintf (fpr,"%g %g %g\n",lp_(0),hp_(1),hp_(2));
    fprintf (fpr,"%g %g %g\n",hp_(0),hp_(1),hp_(2));
  
    if (full) {
      fprintf (fpr,"1 1 1 0\n");
    }
  }

}

