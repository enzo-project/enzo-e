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
  v_[0] = 0;
  v_[1] = 0;
  v_[2] = 0;
}

//----------------------------------------------------------------------

Index Index::index_parent (int min_level) const
{
  Index index = *this;
  int level = index.level();
  if (level > min_level) {
    index.set_child (level,0,0,0);
    index.set_level(level - 1);
  } else {
    ERROR2 ("Index::index_parent",
	    "level %d is less than min_level %d",
	    level,min_level);
  }
  return index;
}

//----------------------------------------------------------------------

Index Index::index_child (int icx, int icy, int icz) const
{
  Index index = *this;
  int level = index.level();
  index.set_level(level+1);
  index.set_child (level+1,icx,icy,icz);
  return index;
}

//----------------------------------------------------------------------

bool Index::is_on_boundary (int axis, int face, int narray) const
{

  int level = this->level();

  // printf ("%s:%d level = %d\n",__FILE__,__LINE__,level);

  ASSERT2 ("Index::is_on_boundary",
	   "%s called for sub-root block level %d",
	   bit_string(5,2,0).c_str(),level,
	   level >= 0);
  
  unsigned array = a_[axis].array;
  unsigned tree  = a_[axis].tree;

  // update tree bits

  unsigned shift_level = (1 << (INDEX_BITS_TREE - level));

  tree += face*shift_level; 

  unsigned shift_overflow = (1 << INDEX_BITS_TREE);

  bool retval = false;

  if (tree & shift_overflow) {

    const int i = array + face;

    retval = ! (0 <= i && i < narray);

  } else {
    
    retval = false;
  }

  return retval;
}

//----------------------------------------------------------------------

bool Index::is_on_boundary 
(
 const int if3[3],
 const int n3[3]
 ) const
{
  return is_on_boundary(0,if3[0],n3[0]) 
    ||   is_on_boundary(1,if3[1],n3[1])
    ||   is_on_boundary(2,if3[2],n3[2]);
}

//----------------------------------------------------------------------

Index Index::index_neighbor (const int if3[3], const int n3[3]) const
{
  Index index = *this;

  const int level = index.level();

  for (int axis = 0; axis < 3; axis++) {
    
    unsigned array = index.a_[axis].array;
    unsigned tree  = index.a_[axis].tree;

    // update tree bits

    unsigned shift_level = (1 << (INDEX_BITS_TREE - level));

    tree += if3[axis]*shift_level; 

    // update array if necessary

    unsigned shift_overflow = (1 << INDEX_BITS_TREE);

    if (tree & shift_overflow) {

      tree &= ~(shift_overflow);

      const int i = array+if3[axis];

      array = (n3[axis] + i) % n3[axis];

    }
    index.a_[axis].array = array;
    index.a_[axis].tree  = tree;
  }
  return index;
}

//----------------------------------------------------------------------

void Index::child 
(int level, int * icx, int * icy, int * icz, int min_level) const
{
  if (level > 0) {
    unsigned shift = INDEX_BITS_TREE - level;
    if (icx) (*icx) = (a_[0].tree >> shift) & 1;
    if (icy) (*icy) = (a_[1].tree >> shift) & 1;
    if (icz) (*icz) = (a_[2].tree >> shift) & 1;
  } else if (level > min_level) {
    unsigned shift = - level;
    if (icx) (*icx) = (a_[0].array >> shift) & 1;
    if (icy) (*icy) = (a_[1].array >> shift) & 1;
    if (icz) (*icz) = (a_[2].array >> shift) & 1;
  }
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
  const unsigned nb = 1 << INDEX_BITS_LEVEL;
  int lx = a_[0].level;
  int ly = a_[1].level;
  int lz = a_[2].level >> 1;
  int sign = (a_[2].level & 1);
  int level = lx + nb*(ly + nb*lz);
  // printf ("%s:%d level lx %d ly %d lz %d sign %d level %d\n",
  // 	  __FILE__,__LINE__,lx,ly,lz,sign,level);
  return sign ? -level : level;
}

//----------------------------------------------------------------------

void Index::set_level(int level)
{ 
  unsigned shift = INDEX_BITS_LEVEL;
  unsigned mask  = ~(1 << shift);
  unsigned sign_mask = level < 0 ? 1 : 0;
  int      sign      = level < 0 ? -1 : 1;
  int      abs_level = level*sign;

  a_[0].level = (abs_level >> (0*shift)) & mask;
  a_[1].level = (abs_level >> (1*shift)) & mask;
  a_[2].level = ((abs_level >> (2*shift-1)) & (mask-1)) | sign_mask;
  clean_();
  ASSERT3 ("Index::set_level",
	   "%s set_level() failed: level %d != level() %d",
	   bit_string(5,3,0).c_str(),level,this->level(),
	   level == this->level());
  
}

//----------------------------------------------------------------------

void Index::set_array(int ix, int iy, int iz)
{ 
  // right-bits = array
  ASSERT4 ("Index::set_array",
	   "Array size (%d %d %d) out of range (maximum %d)",
	   ix,iy,iz,(1<<INDEX_BITS_ARRAY) - 1,
	   ( (ix >> INDEX_BITS_ARRAY) == 0 && 
	     (iy >> INDEX_BITS_ARRAY) == 0 && 
	     (iz >> INDEX_BITS_ARRAY) == 0));

  a_[0].array = ix;
  a_[1].array = iy;
  a_[2].array = iz;
}

//----------------------------------------------------------------------

void Index::tree (int * bx, int * by, int *bz) const
{
  if (bx) (*bx) = a_[0].tree; 
  if (by) (*by) = a_[1].tree; 
  if (bz) (*bz) = a_[2].tree; 
}
  
//----------------------------------------------------------------------

void Index::set_child(int level, int ix, int iy, int iz, int min_level)
{
  if (level > 0) {
    ASSERT ("Index::set_child","level must be at least 1",level>0);
    unsigned shift = (INDEX_BITS_TREE - level);
    unsigned mask  = ~(1 << shift);
    a_[0].tree = (a_[0].tree & mask) | (ix<<shift);
    a_[1].tree = (a_[1].tree & mask) | (iy<<shift);
    a_[2].tree = (a_[2].tree & mask) | (iz<<shift);
  } else if (level > min_level) {
    // switch to array bits if tree bits exhausted (subroot blocks only)
    unsigned shift = - level;
    unsigned mask  = ~(1 << shift);
    a_[0].array = (a_[0].array & mask) | (ix<<shift);
    a_[1].array = (a_[1].array & mask) | (iy<<shift);
    a_[2].array = (a_[2].array & mask) | (iz<<shift);
  }
}

//----------------------------------------------------------------------

void Index::print (const char * msg,
		   int max_level,
		   int rank,
		   bool no_nl,
		   void * simulation ) const
{
  print_(stdout,msg,max_level,rank,no_nl);

#ifdef CELLO_DEBUG
  FILE * fp_debug = ((Simulation *)simulation)->fp_debug();
  print_(fp_debug,msg,max_level,rank,no_nl);
#endif
}

//----------------------------------------------------------------------

void Index::print_ (FILE * fp,
		    const char * msg,
		    int max_level,
		    int rank,
		    bool no_nl) const
{

  std::string buffer;

  const int level = this->level();

  if (max_level == -1) max_level = level;

  int nb = 0;

  for (int axis=0; axis<rank; axis++) {
    nb = std::max(nb,num_bits_(a_[axis].array));
  }

  buffer = buffer + "[ ";
  for (int axis=0; axis<rank; axis++) {

    for (int i=nb; i>=0; i--) {
      int bit = (a_[axis].array & ( 1 << i));
      if (fp != NULL) buffer = buffer + (bit?"1":"0");
    }

    for (int i=0; i<max_level; i++) {
      if (i==0) {
	buffer = buffer + ":";
      }

      if (i < level) {
	int ic3[3];
	child (i+1, &ic3[0], &ic3[1], &ic3[2]);
	buffer = buffer + (ic3[axis] ? "1":"0");
      } else {
	buffer = buffer + " ";
      }
	
    }
    buffer = buffer + " ";
      
  }
  buffer = buffer + "] ";

  buffer = buffer + msg;

  if (! no_nl) buffer = buffer + "\n";

  fprintf (fp,"%s",buffer.c_str());
  fflush(fp);
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
    
  // fprintf (fp,"INDEX %p %s: ", this,msg);
  fprintf (fp,"INDEX %s: ", msg);

  int nb = 0;

  for (int axis=0; axis<rank; axis++) {
    nb = std::max(nb,num_bits_(a_[axis].array));
  }

  for (int axis=0; axis<rank; axis++) {

    for (int i=nb; i>=0; i--) {
      int bit = (a_[axis].array & ( 1 << i));
      fprintf (fp,"%d",bit?1:0);
    }

    for (int level=0; level<max_level; level++) {

      if (level == 0) fprintf (fp,":");
      int ic3[3];
      child (level+1, &ic3[0], &ic3[1], &ic3[2]);
      fprintf (fp,"%d",ic3[axis]);
	
    }
    fprintf (fp," ");
      
  }
  fprintf (fp,"\n");

  fflush(fp);
}

//----------------------------------------------------------------------

std::string Index::bit_string(int max_level,int rank, int num_bits) const
{
  const int level = this->level();

  if (max_level == -1) max_level = this->level();

  std::string bits = "";
  const std::string separator = "_";

  if (num_bits == 0) {
    for (int axis=0; axis<rank; axis++) {
      num_bits = std::max (num_bits,num_bits_(a_[axis].array));
    }
  } else --num_bits;

  for (int axis=0; axis<rank; axis++) {

    for (int i=num_bits; i>=0; i--) {
      int bit = (a_[axis].array & ( 1 << i));
      bits = bits + (bit?"1":"0");
    }

    for (int i=0; i<max_level; i++) {

      if (i == 0) bits = bits + ":";
      if (i < level) {
	int ic3[3];
	child (i+1, &ic3[0], &ic3[1], &ic3[2]);
	bits = bits + (ic3[axis]?"1":"0");
      } else 
	bits = bits + separator;
    }
    if (axis<rank-1) bits = bits + separator;
      
  }
  return bits;
}

//======================================================================

int Index::num_bits_(int value) const
{
  int nb = 32;
  while ( --nb >= 0 && ! (value & (1 << nb))) 
    ;
  return (std::max(nb,0));

}

//----------------------------------------------------------------------

void Index::clean_ ()
{
  for (int level = this->level()+1; level<INDEX_BITS_TREE; level++) {
    set_child(level,0,0,0);
  }
}


