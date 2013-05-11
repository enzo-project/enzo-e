// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Index.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-06
/// @brief    Implementation of Index class for indexing nodes in a forest of octrees

#include "mesh.hpp"

//----------------------------------------------------------------------

Index::Index() { clear(); }

//----------------------------------------------------------------------

Index::Index(const Index & index) 
{
  copy_(index);
}

//----------------------------------------------------------------------

Index::Index(int ix, int iy, int iz) 
{
  clear();
  set_array(ix,iy,iz);
}

//----------------------------------------------------------------------

Index & Index::operator = (const Index & index)
{
  copy_(index);
  return *this;
}

//----------------------------------------------------------------------

bool Index::operator == (const Index & index) const
{
  return (v_[0] == index.v_[0] && 
	  v_[1] == index.v_[1] &&
	  v_[2] == index.v_[2]);
}

//----------------------------------------------------------------------

bool Index::operator != (const Index & index) const
{
  return ! (*this == index);
}

//----------------------------------------------------------------------

// #ifdef CONFIG_USE_CHARM
// #ifndef TEST
// void Index::pup(PUP::er &p) {
//   p(v_,3);
// }
// #endif
// #endif /* CONFIG_USE_CHARM */

// //----------------------------------------------------------------------

void Index::clear () 
{
  for (int axis=0; axis<3; axis++) {
    a_[axis].array = 0;
    a_[axis].tree = 0;
    a_[axis].level = 0;
  }
}

//----------------------------------------------------------------------

Index Index::index_parent () const
{
  TRACE1("index_parent level %d",level());
  Index index = *this;
  int level = index.level();
  if (level > 0) {
    index.set_child (level,0,0,0);
    index.set_level(level - 1);
  } else {
    WARNING("Index::index_parent()",
	    "Attempting to access parent of root");
  }
  return index;
}

//----------------------------------------------------------------------

Index Index::index_child (int ic3[3]) const
{
  TRACE4("index_child level %d  %d %d %d",level(),ic3[0],ic3[1],ic3[2]);
  Index index = *this;
  int level = index.level();
  index.set_level(level+1);
  index.set_child (level+1,ic3[0],ic3[1],ic3[2]);
  return index;
}

//----------------------------------------------------------------------

Index Index::index_neighbor (int axis, int face, int narray) const
{
  TRACE3("index_neighbor axis %d  face %d  narray %d",
	 axis,face,narray);
  if (face == 0) face = -1;

  Index index = *this;

  int array = index.a_[axis].array;
  int tree  = index.a_[axis].tree;

  // update tree bits

  int level = index.level();
  int shift_level = (1 << (INDEX_MAX_TREE_BITS - level));
  TRACE2("shift level %d %0X",level,shift_level);

  tree += face*shift_level; 

  int shift_overflow = (1 << INDEX_MAX_TREE_BITS);

  if (tree & shift_overflow) {

    tree &= ~(shift_overflow);

    TRACE2("array change= %d %d",array,(narray + array + face) % narray);

    array = (narray + array + face) % narray;

  }

  index.a_[axis].array = array;
  index.a_[axis].tree = tree;

  return index;
}

//----------------------------------------------------------------------

Index Index::index_neighbor (int ix, int iy, int iz, int n3[3]) const
{
  TRACE6("index_neighbor ix iy iz  %d %d %d  n3 %d %d %d",
	 ix,iy,iz,n3[0],n3[1],n3[2]);

  Index index = *this;

  int i3[3] = {ix,iy,iz};

  const int level = index.level();

  for (int axis = 0; axis < 3; axis++) {
    
    int array = index.a_[axis].array;
    int tree  = index.a_[axis].tree;

    // update tree bits

    int shift_level = (1 << (INDEX_MAX_TREE_BITS - level));

    tree += i3[axis]*shift_level; 

    // update array if necessary

    int shift_overflow = (1 << INDEX_MAX_TREE_BITS);

    if (tree & shift_overflow) {

      tree &= ~(shift_overflow);

      TRACE2("array change= %d %d",array,(n3[axis] + array + i3[axis]) % n3[axis]);

      array = (n3[axis] + array + i3[axis]) % n3[axis];

    }
    index.a_[axis].array = array;
    index.a_[axis].tree  = tree;
  }


  return index;
}

//----------------------------------------------------------------------

Index Index::index_uncle (int axis, int face, int narray) const
{
  TRACE3("index_uncle axis %d  face %d  narray %d",
	 axis,face,narray);
  return index_parent().index_neighbor(axis,face, narray);
}

//----------------------------------------------------------------------

Index Index::index_nibling (int axis, int face, int ic3[3], int narray) const
{
  TRACE6("index_nibling axis %d  face %d   child %d %d %d  narray %d",
	 axis,face,ic3[0],ic3[1],ic3[2],narray);

  // want facing child in neighbor not corresponding child

  // ic3[axis] = 1 - ic3[axis];

  return index_neighbor(axis,face,narray).index_child(ic3);
}

//----------------------------------------------------------------------

Index Index::index_nibling (int ix, int iy, int iz, int ic3[3], int n3[3]) const
{
  TRACE9("index_nibling ix iy iz %d %d %d  child %d %d %d  n3 %d %d %d",
	 ix,iy,iz,ic3[0],ic3[1],ic3[2],n3[0],n3[1],n3[2]);

  // want facing child in neighbor not corresponding child

  //  ic3[axis] = 1 - ic3[axis];

  return index_neighbor(ix,iy,iz,n3).index_child(ic3);
}

//----------------------------------------------------------------------

void Index::child (int level, int * icx, int * icy, int * icz) const
{
  ASSERT ("Index::child","level must be at least 1",level>0);
  int shift = INDEX_MAX_TREE_BITS - level;
  if (icx) (*icx) = (a_[0].tree >> shift) & 1;
  if (icy) (*icy) = (a_[1].tree >> shift) & 1;
  if (icz) (*icz) = (a_[2].tree >> shift) & 1;
}

//----------------------------------------------------------------------

bool Index::is_root() const
{ return (v_[0]==0 && v_[1]==0 && v_[2]==0); }

//----------------------------------------------------------------------

void Index::array (int * ix, int *iy, int *iz) const
{ 
  if (ix) (*ix) = a_[0].array;
  if (iy) (*iy) = a_[1].array;
  if (iz) (*iz) = a_[2].array;
}

//----------------------------------------------------------------------

int Index::level () const
{
  return a_[0].level + INDEX_MAX_LEVEL_AXIS_RANGE*
    (    a_[1].level + INDEX_MAX_LEVEL_AXIS_RANGE*
	 a_[2].level );  
}

//----------------------------------------------------------------------

void Index::set_level(int level)
{ 
  int shift = INDEX_MAX_LEVEL_AXIS_BITS;
  int mask  = INDEX_MAX_LEVEL_AXIS_RANGE - 1;
  a_[0].level = (level >> (0*shift)) & mask;
  a_[1].level = (level >> (1*shift)) & mask;
  a_[2].level = (level >> (2*shift)) & mask;
  clean();
}

//----------------------------------------------------------------------

void Index::clean ()
{
  for (int level = this->level()+1; level<INDEX_MAX_TREE_BITS; level++) {
    set_child(level,0,0,0);
  }
}

//----------------------------------------------------------------------

void Index::set_array(int ix, int iy, int iz)
{ 
  // right-bits = array
#ifdef CHECK_BOUNDS
  if ( !((0 <= ix && ix < INDEX_MAX_ARRAY_INDEX) ||
	 (0 <= iy && iy < INDEX_MAX_ARRAY_INDEX) ||
	 (0 <= iz && iz < INDEX_MAX_ARRAY_INDEX)) ) {
    printf ("%s:%d  ERROR: out of range (%d %d %d)\n",
	    __FILE__,__LINE__,ix,iy,iz);
    //      return -1;
  }
#endif
  a_[0].array = ix;
  a_[1].array = iy;
  a_[2].array = iz;
}

//----------------------------------------------------------------------

int Index::tree (int axis) const
{ return a_[axis].tree; }
  
//----------------------------------------------------------------------

void Index::set_child(int level, int ix, int iy, int iz)
{
  ASSERT ("Index::set_child","level must be at least 1",level>0);
  int shift = (INDEX_MAX_TREE_BITS - level);
  int mask  = ~(1 << shift);
  a_[0].tree = (a_[0].tree & mask) | (ix<<shift);
  a_[1].tree = (a_[1].tree & mask) | (iy<<shift);
  a_[2].tree = (a_[2].tree & mask) | (iz<<shift);
}

//----------------------------------------------------------------------

void Index::print (const char * msg,
		   int max_level) const
{
  if (max_level == -1) max_level = this->level();

  printf ("INDEX %s: ", msg);

  printf ("A [ ");
  for (int axis=0; axis<3; axis++) {
    printf ("%d ",a_[axis].array);
  }
  printf ("] ");

  printf ("T [ ");
  int ic3[3];
  for (int axis=0; axis<3; axis++) {
    for (int level = 0; level < max_level; level++) {
      child (level+1, &ic3[0], &ic3[1], &ic3[2]);
      printf ("%d",ic3[axis]);
    }
    printf (" ");
  }
  printf ("] ");

  printf ("L [ %d ] ",a_[0].level + INDEX_MAX_LEVEL_AXIS_RANGE*
	  (           a_[1].level + INDEX_MAX_LEVEL_AXIS_RANGE*
		      a_[2].level ));

  printf ("[%08X-%08X-%08X]",v_[0],v_[1],v_[2]);
  printf ("\n");
}

//======================================================================

