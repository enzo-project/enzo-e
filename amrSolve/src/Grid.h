#ifndef GRID_H
#define GRID_H

//======================================================================
//
// File:        Grid.h
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
// History:     2004-03-17  Created
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//======================================================================

#include <stdio.h>

#include "Point.h"

//======================================================================
// Grid class declaration
//======================================================================

class Grid {

  friend class Hierarchy;
  friend class Level;

 public:

  inline Grid () ;
  inline Grid (const Grid &) ;
  inline Grid & operator = (const Grid &) ;
  inline ~Grid () ;

  void attach (class grid *) ;
  void detach () ;

  inline void size (int &, int &, int &) const ;
  inline void range (Point &, Point &) const ;
  inline bool isNeighbor (const Grid &) const ;

  inline bool isLocal () const ;
  inline void print () const ;

  void geomview (FILE *fpr, bool full=true) const ;

 private:

  grid * grid_;
  int n_[3];
  Point lp_, hp_;  // Low-point and high-point

  class Hierarchy * phierarchy_;
  class Level     * plevel_;

};

//======================================================================
// Grid class implementation
//======================================================================

inline Grid::Grid () 
  :  grid_       (0), 
     lp_         (), 
     hp_         (), 
     phierarchy_ (0),
     plevel_     (0)
{ 
  n_[0] = 0; 
  n_[1] = 0; 
  n_[2] = 0; 
}

//----------------------------------------------------------------------

inline Grid::Grid (const Grid &g)               
  :  grid_       (g.grid_) ,
     lp_         (g.lp_) , 
     hp_         (g.hp_), 
     phierarchy_ (g.phierarchy_),
     plevel_     (g.plevel_)
{ 
  n_[0] = g.n_[0];
  n_[1] = g.n_[1];
  n_[2] = g.n_[2];
}

//----------------------------------------------------------------------

inline Grid & Grid::operator = (const Grid &g)  
{ 
  if (this != &g) {
    grid_       = g.grid_;
    n_[0]       = g.n_[0];
    n_[1]       = g.n_[1];
    n_[2]       = g.n_[2];
    lp_         = g.lp_;
    hp_         = g.hp_;
    phierarchy_ = g.phierarchy_;
    plevel_     = g.plevel_;
  }
  return *this;
}

//----------------------------------------------------------------------

inline Grid::~Grid ()
{ 
  grid_ = 0; 
}

//----------------------------------------------------------------------
// Grid.C: void Grid::attach (class grid * grid)
// Grid.C: void Grid::detach ()
//----------------------------------------------------------------------

inline void Grid::size (int & x0, int & x1, int & x2) const 
{ 
  x0 = n_[0];
  x1 = n_[1];
  x2 = n_[2];
}

//----------------------------------------------------------------------

inline void Grid::range (Point & p1, Point & p2) const 
{
  p1 = lp_;
  p2 = hp_;
}

//----------------------------------------------------------------------

inline bool Grid::isNeighbor (const Grid &g) const 
// Returns true iff grid--including one ghost zone layer--overlap
{
  // Grids are not self-neighbors!

  if (this == &g) return false;

  assert (this->plevel_ == g.plevel_);

  // Determine tolerance for comparing positions:
  // use half a grid zone size

  double hh[3];
  hh[0] = 0.5 * (hp_(0)-lp_(0))/n_[0];    
  hh[1] = 0.5 * (hp_(1)-lp_(1))/n_[1];
  hh[2] = 0.5 * (hp_(2)-lp_(2))/n_[2];

  // Grids are adjacent iff they are not far
  // Grids are far iff they are far in at least one axis
  // Grids p and q are far along an axis iff ph < ql or qh < pl
  //   ( lppph lqqqh or lqqqh lppph )
  // To prevent possible float noise, assume p < q iff p < q-hh

  return !
    ((hp_(0) < g.lp_(0) - hh[0]) || (g.hp_(0) < lp_(0) - hh[0])) ||
    ((hp_(1) < g.lp_(1) - hh[1]) || (g.hp_(1) < lp_(1) - hh[1])) ||
    ((hp_(2) < g.lp_(2) - hh[2]) || (g.hp_(2) < lp_(2) - hh[2]));

  
}

//----------------------------------------------------------------------

inline bool Grid::isLocal () const 
{
  fprintf (stderr, "WARNING: Grid::isLocal () not implemented!\n");
  return 0;
}

//----------------------------------------------------------------------

inline void Grid::print () const 
{
  printf ("amrSolve::Grid-ID    %p\n",this);
  printf ("amrSolve::Grid-size  %d %d %d\n", n_[0], n_[1], n_[2]);
  printf ("amrSolve::Grid-range %f %f %f  %f %f %f\n", 
	  lp_(0), lp_(1), lp_(2), hp_(0), hp_(1), hp_(2));
}

#endif
