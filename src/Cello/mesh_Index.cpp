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

Index Index::index_child (int icx, int icy, int icz) const
{
  TRACE4("index_child level %d  %d %d %d",level(),icx,icy,icz);
  Index index = *this;
  int level = index.level();
  index.set_level(level+1);
  index.set_child (level+1,icx,icy,icz);
#ifdef CELLO_TRACE
  char buffer[40];
  sprintf (buffer,"index_child(ix iy iz %d %d %d)",icx,icy,icz);
  this->print(buffer);
  index.print(buffer);
#endif
  return index;
}

//----------------------------------------------------------------------

Index Index::index_neighbor (int axis, int face, int narray, bool periodic) const
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

      // if (!periodic) {
      //   // return self index if on boundary
      //   if ( ! ((0 <= array+face) && (array+face < narray))) return index;
      // }

    array = (narray + array + face) % narray;

    TRACE3("array %d  face %d  narray %d",array,face,narray);
  }

  index.a_[axis].array = array;
  index.a_[axis].tree = tree;

#ifdef CELLO_TRACE
  char buffer[40];
  sprintf (buffer,"index_neighbor(axis %d face %d)",axis,face);
  this->print(buffer);
  index.print(buffer);
#endif

  return index;
}

//----------------------------------------------------------------------

bool Index::is_on_boundary (int axis, int face, int narray, bool periodic) const
{
  if (periodic) return false;

  TRACE4("is_on_boundary axis %d  face %d  narray %d period %d", 
	 axis,face,narray,periodic);

  if (face == 0) face = -1;

  int level = this->level();
  int array = a_[axis].array;
  int tree  = a_[axis].tree;

  // update tree bits

  int shift_level = (1 << (INDEX_MAX_TREE_BITS - level));
  TRACE2("shift level %d %0X",level,shift_level);

  tree += face*shift_level; 

  int shift_overflow = (1 << INDEX_MAX_TREE_BITS);

  if (tree & shift_overflow) {

    if (narray == 1) return true;

    tree &= ~(shift_overflow);

    TRACE2("array change= %d %d",array,(narray + array + face) % narray);

    return ! ((0 <= array+face) && (array+face < narray) );

  } else {
    
    return false;
  }

}

//----------------------------------------------------------------------

Index Index::index_neighbor (int if3[3], int n3[3], bool periodic) const
{
  return index_neighbor(if3[0],if3[1],if3[2],n3,periodic);
}
//----------------------------------------------------------------------

Index Index::index_neighbor (int ix, int iy, int iz, int n3[3], bool periodic) const
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

#ifdef CELLO_TRACE
  char buffer[40];
  sprintf (buffer,"index_neighbor(ix iy iz %d %d %d)",ix,iy,iz);
  this->print(buffer);
  index.print(buffer);
#endif

  return index;
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
		   int max_level,
		   int rank) const
{
  if (max_level == -1) max_level = this->level();
    
  PARALLEL_PRINTF ("INDEX %p %s: ", this,msg);

  int nb = 0;

  for (int axis=0; axis<rank; axis++) {
    nb = std::max(nb,num_bits_(a_[axis].array));
  }

  for (int axis=0; axis<rank; axis++) {

    print_bits_(a_[axis].array,nb);

    PARALLEL_PRINTF (":");

    for (int level=0; level<max_level; level++) {

      int ic3[3];
      child (level+1, &ic3[0], &ic3[1], &ic3[2]);
      PARALLEL_PRINTF ("%d",ic3[axis]);
	
    }
    PARALLEL_PRINTF (" ");
      
  }
  PARALLEL_PRINTF ("\n");

  fflush(stdout);
}

//----------------------------------------------------------------------

void Index::write (int ip,
		   const char * msg,
		   int max_level,
		   int rank) const
{
  char filename[80];
  sprintf (filename,"index.%s.%d",msg,ip);
  FILE * fp = fopen(filename,"a");

  if (max_level == -1) max_level = this->level();
    
  fprintf (fp,"INDEX %p %s: ", this,msg);

  int nb = 0;

  for (int axis=0; axis<rank; axis++) {
    nb = std::max(nb,num_bits_(a_[axis].array));
  }

  for (int axis=0; axis<rank; axis++) {

    write_bits_(fp,a_[axis].array,nb);

    fprintf (fp,":");

    for (int level=0; level<max_level; level++) {

      int ic3[3];
      child (level+1, &ic3[0], &ic3[1], &ic3[2]);
      fprintf (fp,"%d",ic3[axis]);
	
    }
    fprintf (fp," ");
      
  }
  fprintf (fp,"\n");

  fflush(stdout);
}

//----------------------------------------------------------------------

int Index::num_bits_(int value) const
{
  int nb = 32;
  while ( --nb >= 0 && ! (value & (1 << nb))) 
    ;
  return (std::max(nb,0));

}

//----------------------------------------------------------------------

void Index::print_bits_(int value, int nb) const
{
  for (int i=nb; i>=0; i--) {
    int bit = (value & ( 1 << i));
    PARALLEL_PRINTF ("%d",bit?1:0);
  }
}

//----------------------------------------------------------------------

void Index::write_bits_(FILE * fp,int value, int nb) const
{
  for (int i=nb; i>=0; i--) {
    int bit = (value & ( 1 << i));
    fprintf (fp,"%d",bit?1:0);
  }
}

//----------------------------------------------------------------------

std::string Index::bit_string(int max_level,int rank) const
{
  std::string bits = "";

  int nba = 0;
  for (int axis=0; axis<rank; axis++) {
    nba = std::max(nba,num_bits_(a_[axis].array));
  }

  for (int axis=0; axis<rank; axis++) {

    for (int i=nba; i>=0; i--) {
      int bit = (a_[axis].array & ( 1 << i));
      bits = bits + (bit?"1":"0");
    }

    bits = bits + ":";

    for (int level=0; level<max_level; level++) {

      int ic3[3];
      child (level+1, &ic3[0], &ic3[1], &ic3[2]);
      bits = bits + (ic3[axis]?"1":"0");
	
    }
    if (axis<rank-1) bits = bits + " ";
      
  }
  return bits;
}

//======================================================================

