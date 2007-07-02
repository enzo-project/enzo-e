#ifndef AMR_ITLEVELS_H
#define AMR_ITLEVELS_H

//======================================================================
//
// File:        ItLevels.h
//
// Package:     amrSolve
//
// Description: Class for iterating through Levels in a Hierarchy
//
// Classes:     ItLevels
//
//----------------------------------------------------------------------
//
// Author:      James Bordner (jbordner@cosmos.ucsd.edu)
//
// History:     2004-03-24  Created
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//======================================================================

//======================================================================
// ItLevels class declaration
//======================================================================
  
class ItLevels {
    
 public:
    
  inline ItLevels (const Hierarchy &);
  inline ItLevels (const ItLevels &);
  inline ItLevels & operator = (const ItLevels &);
  inline ~ItLevels ();
    
  inline const Level * operator ++ ();
  inline const Level * operator *() const;
  inline const Level * begin () const;
  inline const Level * end () const;
  inline void reset () ;
    
 private:
    
  const Hierarchy * hierarchy_;
  int ilevel_;  // 0 <= ilevel_ <= hierarchy_.numLevels();
    
};
  
//======================================================================
// ItLevels class implementation
//======================================================================

inline ItLevels::ItLevels (const Hierarchy &H)
  : hierarchy_ (&H),
    ilevel_ (H.numLevels())
{
}

//----------------------------------------------------------------------

inline ItLevels::ItLevels (const ItLevels &it)
  : hierarchy_ (it.hierarchy_),
    ilevel_ (it.ilevel_)
{
}

//----------------------------------------------------------------------

inline ItLevels & ItLevels::operator = (const ItLevels &it)
{
  if (this != &it) {
    hierarchy_ = it.hierarchy_;
    ilevel_ = it.ilevel_;
  }
  return *this;
}

//----------------------------------------------------------------------

inline ItLevels::~ItLevels ()
{
  hierarchy_ = 0;
  ilevel_ = 0;
}

//----------------------------------------------------------------------
    
inline const Level * ItLevels::operator ++ ()
{
  ilevel_ = (ilevel_ + 1) % (hierarchy_->numLevels() + 1) ;
  return *(*this) ;
}

//----------------------------------------------------------------------

inline const Level * ItLevels::operator *() const
{
  return (ilevel_ < hierarchy_->numLevels()) ? hierarchy_->level(ilevel_) : 0;
}

//----------------------------------------------------------------------

inline const Level * ItLevels::begin () const
{
  return hierarchy_->level(0);
}

//----------------------------------------------------------------------

inline const Level * ItLevels::end () const
{
  return 0;
}

//----------------------------------------------------------------------

inline void ItLevels::reset () 
{
  ilevel_ = hierarchy_->numLevels();
}

//----------------------------------------------------------------------

#endif
