#ifndef INDEX_HPP
#define INDEX_HPP

#include "error_Error.hpp"

// ARRAY_BITS + TREE_BITS + LEVEL_AXIS_BITS == 32

#define INDEX_MAX_ARRAY_BITS       10
#define INDEX_MAX_ARRAY_INDEX    1024 /* 2 ** INDEX_MAX_ARRAY_BITS */
#define INDEX_MAX_TREE_BITS        20
#define INDEX_MAX_LEVEL_AXIS_BITS   2
#define INDEX_MAX_LEVEL_AXIS_RANGE  4 /* 2 ** AXIS_BITS2 */

struct BIndex {

  // maximum 1024 x 1024 x 1024 root nodes

  unsigned array : INDEX_MAX_ARRAY_BITS;

  // maximum INDEX_MAX_TREE_BITS levels / bits
  // L    T
  // 0    NULL
  // 1    0--  | 0-- | 0-- | 0-- | 1-- | 1-- | 1-- | 1-- |
  // 2    00-  | 00- | 01- | 01- | 10- | 10- | 11- | 11- |
  // 3    000  | 001 | 010 | 011 | 100 | 101 | 110 | 111 |
  //  bit fields left-justified to simplify finding neighbors
  //  - : unaccessed, but must be 0
  //  INDEX_MAX_TREE_BITS - LEVEL << 1 is lowest bit for LEVEL

  unsigned  tree : INDEX_MAX_TREE_BITS; 

  // maximum 32 level specification

  unsigned level : INDEX_MAX_LEVEL_AXIS_BITS; 

};

class Index {
  //
  // three integers a_[0] a_[1] a_[2] one for each axis
  // 96 bits total
  // Array index 10x3   30 1024
  // Tree  index 20x3   60 20 levels
  // Level index  6x1   32      
  //       [       |       |       |        )
  // a_[0] [TTTTTTTTTTTTTTTTTTTTLLAAAAAAAAAA)
  // a_[1] [TTTTTTTTTTTTTTTTTTTTLLAAAAAAAAAA)
  // a_[2] [TTTTTTTTTTTTTTTTTTTTLLAAAAAAAAAA)
  // ABITS  10  one index per integer
  // TBITS  20  L bit indices per integer
  // LBITS   2  6 bits spread across three integers
public:

  Index();

  Index(const Index & index);

  Index(int ix, int iy, int iz);

  Index & operator = (const Index & index);

  bool operator == (const Index & index) const;

  bool operator != (const Index & index) const;

  void clear () ;
  
  Index index_parent () const;

  Index index_child (int ic3[3]) const
  { return index_child(ic3[0],ic3[1],ic3[2]); }

  Index index_child (int icx, int icy, int icz) const;

  /// Index of neighbor along given face
  Index index_neighbor 
  (int axis, int face, int narray, bool periodic=false) const;

  /// Generalized neighbor including edge and corner
  Index index_neighbor 
  (int ix, int iy, int iz, int n3[3], bool periodic=false) const;

  inline Index index_uncle 
  (int axis, int face, int narray, bool periodic=false) const
  {  return index_parent().index_neighbor(axis,face, narray, periodic); }

  /// return the given parent's neighbor index
  inline Index index_nibling 
  (int axis, int face, int ic3[3], int narray,bool periodic=false) const
  {  return index_neighbor(axis,face,narray,periodic).index_child(ic3); }

  /// Return the given neighbor's child index
  Index index_nibling 
  (int ix, int iy, int iz, int ic3[3], int n3[3],bool periodic=false) const
  {  return index_neighbor(ix,iy,iz,n3,periodic).index_child(ic3); }

  /// Whether the face is on the domain boundary
  bool is_on_boundary 
  (int axis, int face, int narray, bool periodic=false) const;

  /// child index of this node in parent
  void child (int level, int * icx, int * icy, int * icz) const;

  /// Return whether this is the "root" node in the forest
  /// (array (0 0 0), level 0)
  bool is_root() const;

  /// Return the indices of the level-0 node containing this node
  void array (int * ix, int *iy, int *iz) const;

  /// Return the level of this node
  int level () const;

  /// Return the packed bit index for the given axis
  // unsigned value (int axis) const;

  /// Set the Index according to raw bit values
  void set_values (int v3[3])
  {
    v_[0] = v3[0];
    v_[1] = v3[1];
    v_[2] = v3[2];
  }

  /// Return the packed bit index for the given axis
  void values (int v3[3]) const
  { v3[0] = v_[0];
    v3[1] = v_[1];
    v3[2] = v_[2];
  }


  /// Set the level for this node
  void set_level(int level);

  /// Clear tree bits that are associated with levels higher than
  /// the actual level
  void clean ();

  /// Accumulate array part of an index
  void set_array(int ix, int iy, int iz);

  /// Return the packed tree bits for the given axis
  int tree (int axis) const;
  
  /// Set the child indicies of this node in the parent
  void set_child(int level, int ix, int iy=0, int iz=0);

  void print (const char * msg = "\0",
	      int max_level = -1,
	      int rank = 3) const;

  std::string bit_string (int max_level,int rank) const;

private: // functions

  int num_bits_(int value) const;
  void print_bits_(int value, int nb) const;

  inline void copy_ (const Index & index)
  {
    a_[0] = index.a_[0];
    a_[1] = index.a_[1];
    a_[2] = index.a_[2];
  }

private: // attributes

    union {
      BIndex a_[3];
      int v_[3];
    };
};

#ifdef CONFIG_USE_CHARM
#ifndef TEST
  PUPbytes(Index);
#endif
#endif /* CONFIG_USE_CHARM */

#ifdef CONFIG_USE_CHARM
#ifndef TEST
// public:
//   void pup(PUP::er &p) {
    
//   }
PUPbytes(BIndex);
#endif
#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------
#ifdef CONFIG_USE_CHARM

#ifndef TEST
class CkArrayIndexIndex:public CkArrayIndex {
  Index index_;
public:
  CkArrayIndexIndex(const Index &in)
  {
    index_=in;
    nInts=sizeof(index_)/sizeof(int);
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    // ADDED SINCE OTHERWISE CkArrayIndex::index[] DOES NOT GET INITIALIZED
    int v3[3];
    in.values(v3);
    index[0] = v3[0];
    index[1] = v3[1];
    index[2] = v3[2];
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  }
  //Not required, but convenient: cast-to-foo operators
   operator Index &() {return index_;}
   operator const Index &() const {return index_;}
    Index & ind() { return index_; }
};
#endif
#endif /* INDEX_HPP */

#endif /* CONFIG_USE_CHARM */
