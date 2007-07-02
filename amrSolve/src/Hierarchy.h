#ifndef HIERARCHY_H
#define HIERARCHY_H

//======================================================================
//
// File:        Hierarchy.h
//
// Package:     amrSolve
//
// Description: Interface to Enzo hierarchy
//
// Classes:     Hierarchy
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

#include <vector>
#include <map>

//======================================================================
// Hierarchy class declaration
//======================================================================

class Hierarchy {

 public:

  inline Hierarchy () ;
  inline Hierarchy (const Hierarchy &) ;
  inline Hierarchy & operator = (const Hierarchy &) ;
  inline ~Hierarchy () ;

  void attach (LevelHierarchyEntry *[]) ;
  void detach () ;

  inline int numLevels () const ;
  inline Level * level (int) const ;

  inline void assertParent (const Grid &p, const Grid &c) ;

  void geomview (FILE *fpr, bool full=true) const ;

 private:

  // Array of levels in the hierarchy

  std::vector   <Level *>        levels_ ;

  // Each grid's parent

  std::map      <const Grid *, const Grid *> parents_;

  // Each grid's list of children

  std::multimap <const Grid *, const Grid *> children_;

};

//======================================================================
// Hierarchy class implementation
//======================================================================

//----------------------------------------------------------------------

inline Hierarchy::Hierarchy () 
  : levels_()
{
}

//----------------------------------------------------------------------

inline Hierarchy::Hierarchy (const Hierarchy &hierarchy)              
{
  levels_=hierarchy.levels_ ;
}

//----------------------------------------------------------------------

inline Hierarchy & Hierarchy::operator = (const Hierarchy &hierarchy) 
{
  if (this != &hierarchy) {
    levels_ = hierarchy.levels_ ; 
  }
  return *this ;
}

//----------------------------------------------------------------------

inline Hierarchy::~Hierarchy ()                        
{
}

//----------------------------------------------------------------------
// Hierarchy.C: void Hierarchy::attach (LevelHierarchyEntry * [])
// Hierarchy.C: void Hierarchy::detach ()
//----------------------------------------------------------------------

inline int Hierarchy::numLevels () const
{ 
  return levels_.size() ; 
}

//----------------------------------------------------------------------

inline Level * Hierarchy::level (int i) const
{
  if (0 <= i && i < levels_.size()) {
    return levels_[i];
  } else {
    return 0;
  }
}

//----------------------------------------------------------------------

inline void Hierarchy::assertParent (const Grid &p, const Grid &c)
{
  parents_.insert (std::make_pair(&c,&p));
  children_.insert(std::make_pair(&p,&c));
}

#endif
