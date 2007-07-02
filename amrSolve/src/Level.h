#ifndef LEVEL_H
#define LEVEL_H

//======================================================================
//
// File:        Level.h
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
// History:     2004-03-17  Created
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//======================================================================

#include "Grid.h"

//======================================================================
// Level class declaration
//======================================================================

class Level {

 public:

  inline Level () ;
  inline Level (const Level &) ;
  inline Level & operator = (const Level &) ;
  inline ~Level () ;

  void attach (LevelHierarchyEntry *) ;
  void detach () ;

  inline int numGrids () const;
  inline Grid * grid (int) const;

  inline void assertNeighbors (const Grid &, const Grid &);

  void geomview (FILE *fpr, bool full=true) const ;

 private:

  // Array of grids in the level

  std::vector   <Grid *>         grids_;     

  // Each grid's list of neighbors
  
  std::multimap <const Grid *, const Grid *> neighbors_; 

};

//======================================================================
// Level class implementation 
//======================================================================

inline Level::Level () 
  : grids_()
{
}

//----------------------------------------------------------------------

inline Level::Level (const Level &level)
  : grids_ (level.grids_)

{
}

//----------------------------------------------------------------------

inline Level & Level::operator = (const Level &level)
{
  if (this != &level) {
    grids_ = level.grids_; 
  }
  return *this;
}

//----------------------------------------------------------------------

inline Level::~Level ()                        
{
}

//----------------------------------------------------------------------
// Level.C: void Level::attach (LevelHierarchyEntry * )
// Level.C: void Level::detach ()
//----------------------------------------------------------------------

inline int Level::numGrids () const
{ 
  return grids_.size(); 
}

//----------------------------------------------------------------------

inline Grid * Level::grid (int i) const
{
  if (0 <= i && i < grids_.size()) {
    return grids_[i];
  } else {
    return 0;
  }
}

//----------------------------------------------------------------------

inline void Level::assertNeighbors (const Grid &g1, const Grid &g2)
{
  if (&g1 != &g2) {
    neighbors_.insert(std::make_pair(&g1,&g2));
    neighbors_.insert(std::make_pair(&g2,&g1));
  }
}

//----------------------------------------------------------------------

#endif
