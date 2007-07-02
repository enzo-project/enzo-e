#ifndef AMR_ITGRIDS_H
#define AMR_ITGRIDS_H

//======================================================================
//
// File:        ItGrids.h
//
// Package:     amrSolve
//
// Description: Class for iterating through Grids in a Level
//
// Classes:     ItGrids
//
//----------------------------------------------------------------------
//
// Author:      James Bordner (jbordner@cosmos.ucsd.edu)
//
// History:     2004-03-23  Created
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//======================================================================

//======================================================================
// ItGrids class declaration
//======================================================================
  
class ItGrids {
    
 public:
    
  inline ItGrids (const Level &);
  inline ItGrids (const ItGrids &);
  inline ItGrids & operator = (const ItGrids &);
  inline ~ItGrids ();
    
  inline const Grid * operator ++ ();
  inline const Grid * operator *() const;
  inline const Grid * begin () const;
  inline const Grid * end () const;
  inline void reset () ;
    
 private:
    
  const Level * level_;
  int igrid_;  // 0 <= igrid_ <= level_.numGrids();
    
};
  
//======================================================================
// ItGrids class implementation
//======================================================================

inline ItGrids::ItGrids (const Level &L)
  : level_ (&L),
    igrid_ (L.numGrids())
{
}

//----------------------------------------------------------------------

inline ItGrids::ItGrids (const ItGrids &it)
  : level_ (it.level_),
    igrid_ (it.igrid_)
{
}

//----------------------------------------------------------------------

inline ItGrids & ItGrids::operator = (const ItGrids &it)
{
  if (this != &it) {
    level_ = it.level_;
    igrid_ = it.igrid_;
  }
  return *this;
}

//----------------------------------------------------------------------

inline ItGrids::~ItGrids ()
{
  level_ = 0;
  igrid_ = 0;
}

//----------------------------------------------------------------------
    
inline const Grid * ItGrids::operator ++ ()
{
  igrid_ = (igrid_ + 1) % (level_->numGrids() + 1) ;
  return *(*this) ;
}

//----------------------------------------------------------------------

inline const Grid * ItGrids::operator *() const
{
  return (igrid_ < level_->numGrids()) ? level_->grid(igrid_) : 0;
}

//----------------------------------------------------------------------

inline const Grid * ItGrids::begin () const
{
  return level_->grid(0);
}

//----------------------------------------------------------------------

inline const Grid * ItGrids::end () const
{
  return 0;
}

//----------------------------------------------------------------------

inline void ItGrids::reset () 
{
  igrid_ = level_->numGrids();
}

//----------------------------------------------------------------------

#endif
