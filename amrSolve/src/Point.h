#ifndef POINT_H
#define POINT_H

//======================================================================
//
// File:        Point.h
//
// Package:     amrSolve
//
// Description: A point in 3D
//
// Classes:     Point
//
//----------------------------------------------------------------------
//
// Author:      James Bordner (jbordner@cosmos.ucsd.edu)
//
// History:     2002-03-18  Created
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//======================================================================

#define Scalar double

//======================================================================
// Point class interface
//======================================================================

class Point {

public:

  inline Point ();
  inline Point (Scalar, Scalar, Scalar);
  inline Point (const Point &);
  inline Point & operator = (const Point &);
  inline ~Point ();

  inline void get (Scalar &, Scalar &, Scalar &) const;
  inline void set (Scalar,   Scalar,   Scalar);

  inline Scalar operator () (int) const;
  inline Scalar & operator () (int);

  friend bool operator <= (const Point &, const Point &);

  friend void sort (Point &, Point &);

private:

  Scalar a_[3];

};

//======================================================================
// Point class implementation
//======================================================================

inline Point::Point () 
{ 
  a_[0] = 0;
  a_[1] = 0;
  a_[2] = 0; 
}

//----------------------------------------------------------------------

inline Point::Point (Scalar a0, Scalar a1, Scalar a2)
{ 
  a_[0] = a0;
  a_[1] = a1;
  a_[2] = a2; 
}

//----------------------------------------------------------------------

inline Point::Point (const Point &p)
{ 
  a_[0] = p.a_[0];
  a_[1] = p.a_[1];
  a_[2] = p.a_[2]; 
}

//----------------------------------------------------------------------

inline Point & Point::operator = (const Point &p)
{ 
  a_[0] = p.a_[0];
  a_[1] = p.a_[1];
  a_[2] = p.a_[2]; 
  return *this; 
}

//----------------------------------------------------------------------

inline Point::~Point ()
{
}

//----------------------------------------------------------------------

inline void Point::get (Scalar &a0, Scalar &a1, Scalar &a2) const
{ 
  a0 = a_[0]; 
  a1 = a_[1]; 
  a2 = a_[2]; 
}

//----------------------------------------------------------------------

inline void Point::set ( Scalar a0,  Scalar a1,  Scalar a2)
{ 
  a_[0] = a0; 
  a_[1] = a1; 
  a_[2] = a2; 
}

//----------------------------------------------------------------------

inline Scalar Point::operator () (int i) const
{ 
  return a_[i]; 
}

//----------------------------------------------------------------------

inline Scalar & Point::operator () (int i)
{ 
  return a_[i]; 
}

//----------------------------------------------------------------------

inline bool operator <= (const Point &p1, const Point &p2)
{ 
  return p1.a_[0] <= p2.a_[0]
    &&   p1.a_[1] <= p2.a_[1]
    &&   p1.a_[2] <= p2.a_[2];
}

//----------------------------------------------------------------------

#define SWAP(a,b) { Scalar t; t = a; a = b; b = t; }

inline void sort (Point &p1, Point &p2)

{ 
  Scalar t;
  if (p1.a_[0] > p2.a_[0]) SWAP(p1.a_[0],p2.a_[0]);
  if (p1.a_[1] > p2.a_[1]) SWAP(p1.a_[1],p2.a_[1]);
  if (p1.a_[2] > p2.a_[2]) SWAP(p1.a_[2],p2.a_[2]); 
}

#endif
