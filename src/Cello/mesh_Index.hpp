#ifndef INDEX_HPP
#define INDEX_HPP

#include "error_Error.hpp"

// BITS_ARRAY + BITS_TREE + BITS_LEVEL == 32

#define INDEX_BITS_ARRAY  10
#define INDEX_BITS_TREE   20
#define INDEX_BITS_LEVEL   2

struct BIndex {

  // maximum 1024 x 1024 x 1024 root nodes

  unsigned array : INDEX_BITS_ARRAY;
  unsigned  tree : INDEX_BITS_TREE; 
  unsigned level : INDEX_BITS_LEVEL; 

  // maximum INDEX_BITS_TREE levels / bits
  // L    T
  // 0    NULL
  // 1    0--  | 0-- | 0-- | 0-- | 1-- | 1-- | 1-- | 1-- |
  // 2    00-  | 00- | 01- | 01- | 10- | 10- | 11- | 11- |
  // 3    000  | 001 | 010 | 011 | 100 | 101 | 110 | 111 |
  //  bit fields left-justified to simplify finding neighbors
  //  - : unaccessed, but must be 0
  //  INDEX_BITS_TREE - LEVEL << 1 is lowest bit for LEVEL
  // level range (-2^5, 2^5) = [-31, 31]
  // level sign bit is v_[2] & (1 << (INDEX_BITS_LEVEL - 1))

};

class Index {
  //
  // three integers a_[0] a_[1] a_[2] one for each axis
  // 96 bits total
  // Array index 10x3   30 1024
  // Tree  index 20x3   60 20 levels
  // Level index  6x1   32 + sign for sub-root blocks      
  //       [       |       |       |        )
  // a_[0] [AAAAAAAAAATTTTTTTTTTTTTTTTTTTTLL)
  // a_[1] [AAAAAAAAAATTTTTTTTTTTTTTTTTTTTLL)
  // a_[2] [AAAAAAAAAATTTTTTTTTTTTTTTTTTTTLS)
  // ABITS  10  one index per integer a_[i]
  // TBITS  20  bit indices per integer a_[i]
  // LBITS   2  6 bits spread across three integers
  //            with one bit representing sign ('S')
public:

  Index();

  Index(const Index & index);

  Index(int ix, int iy, int iz);

  Index & operator = (const Index & index);

  bool operator == (const Index & index) const;

  bool operator != (const Index & index) const;

  void clear () ;
  
  Index index_parent (int min_level = 0) const;

  Index index_child (const int ic3[3]) const
  { return index_child(ic3[0],ic3[1],ic3[2]); }

  Index index_child (int icx, int icy, int icz) const;

  /// Return the index for the given neighbor
  Index index_neighbor (const int if3[3], const int n3[3]) const;

  /// Whether the face is on the domain boundary
  bool is_on_boundary 
  (int axis, int face, int narray) const;

  /// Whether the face is on the domain boundary
  bool is_on_boundary (const int if3[3], const int n3[3]) const;

  /// Return whether this is the "root" node in the forest
  /// (array (0 0 0), level 0)
  bool is_root() const;

  /// Return the level of this node
  int level() const;

  /// Return the packed bit index for the given axis
  // unsigned value (int axis) const;

  /// Set the Index according to raw bit values
  void set_values (const int v3[3])
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

  /// Return the indices of the level-0 node containing this node
  void array (int * ix, int *iy, int *iz) const;

  /// Accumulate array part of an index
  void set_array(int ix, int iy, int iz);

  /// Return the packed tree bits for each axis
  void tree (int * bx = 0, int *by = 0, int *bz = 0) const;
  
  /// child index of this node in parent
  void child (int level, int * ix, int * iy, int * iz, int min_level = 0) const;

  /// Set the child indicies of this node in the parent
  void set_child(int level, int ix, int iy=0, int iz=0, int min_level = 0);

  void print (const char * msg,
	      int max_level,
	      int rank,
	      bool no_nl,
	      void * simulation = 0) const;
  

  void write (int ip,
	      const char * msg = "\0",
	      int max_level = -1,
	      int rank = 3) const;

  std::string bit_string (int max_level,int rank, const int nb3[3]) const;

  /// Comparison operator required for Charm++ pup()
  friend bool operator < (const Index & x, const Index & y) {
    if (x.v_[2] < y.v_[2]) return true;
    if (x.v_[2] > y.v_[2]) return false;
    if (x.v_[1] < y.v_[1]) return true;
    if (x.v_[1] > y.v_[1]) return false;
    return  (x.v_[0] < y.v_[0]);
  }

private: // functions

  /// Clear tree bits that are associated with levels higher than
  /// the actual level
  void clean_ ();

  int num_bits_(int value) const;
	
  inline void copy_ (const Index & index)
  {
    v_[0] = index.v_[0];
    v_[1] = index.v_[1];
    v_[2] = index.v_[2];
  }

  void print_ (FILE * fp,
	       const char * msg,
	       int max_level,
	       int rank,
	       bool no_nl) const;

private: // attributes

    union {
      BIndex a_[3];
      int v_[3];
    };
};

#ifndef TEST
  PUPbytes(Index);
#endif

#ifndef TEST
// public:
//   void pup(PUP::er &p) {
    
//   }
PUPbytes(BIndex);
#endif

//----------------------------------------------------------------------
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
