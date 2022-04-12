#ifndef INDEX_HPP
#define INDEX_HPP

#include "error_Error.hpp"

// BITS_ARRAY + BITS_TREE + BITS_LEVEL == 32

#define INDEX_BITS_ARRAY  10
#define INDEX_BITS_TREE   20
#define INDEX_BITS_LEVEL   2

#define INDEX_UNDEFINED_LEVEL -999
class NodeBits {

  // original order ATL crashed in Charm++ during load balancing

public:

  unsigned level : INDEX_BITS_LEVEL;
  unsigned array : INDEX_BITS_ARRAY;
  unsigned  tree : INDEX_BITS_TREE;

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

  Index(int iax, int iay, int iaz);

  bool operator == (const Index & index) const;

  bool operator != (const Index & index) const;

  inline int operator [] (std::size_t i) const
  { return v_[i]; }

  void clear () ;

  Index index_parent (int min_level = 0) const;

  inline Index index_child (const int ic3[3], int min_level=0) const
  { return index_child(ic3[0],ic3[1],ic3[2],min_level); }

  Index index_child (int icx, int icy, int icz, int min_level = 0) const;

  /// Return the index for the given neighbor
  Index index_neighbor (const int if3[3], const int na3[3]) const;

  /// Return the index of the ancestor in the given level_ancestor <= level
  /// default is root level
  Index index_ancestor (int level_ancestor = 0, int min_level = 0) const;

  /// Return offset of block in given level
  void index_level (int i3[3], int level) const;
  int index_level (int level,int axis) const;

  /// Whether the face is on the domain boundary
  bool is_on_boundary
  (int axis, int face, int narray) const;

  /// Whether the face is on the domain boundary
  bool is_on_boundary (const int if3[3], const int na3[3]) const;

  /// Whether an index is in the same subtree relative to a given
  /// root level
  bool is_in_same_subtree (Index index, int min_level = 0, int root_level = 0);

  /// Return whether this is the "root" node in the array of octrees
  /// (array (0 0 0), level 0)
  bool is_root() const;

  /// Whether given `index` is a sibling (has the same parent as this Index)
  bool is_sibling (Index index) const
  {
    const int level = this->level();
    return (level >= 1 && index.level() >= 1) ?
      (index_parent() == index.index_parent()) : false;
  }
  /// Whether given `index` is a "nibling" (child of a sibling of this Index)
  bool is_nibling (Index index) const
  {
    const int level = this->level();
    return (level >= 1 && index.level() >= 2) ?
      (index_parent() == index.index_parent().index_parent()) : false;
  }

  /// Return the dimensionality of shared face (0 corner, 1 edge, 2
  /// plane), or -1 if disjoint
  int adjacency (Index index, int rank, const int p3[3]) const;

  /// Return refinement level of the Index
  int level() const;

  /// Return the packed bit index for the given axis
  // unsigned value (int axis) const;

  /// Set the Index according to raw bit values
  inline void set_values (const int v3[3])
  {
    v_[0] = v3[0];
    v_[1] = v3[1];
    v_[2] = v3[2];
  }

  /// Return the packed bit index for the given axis
  inline void values (int v3[3]) const
  { v3[0] = v_[0];
    v3[1] = v_[1];
    v3[2] = v_[2];
  }


  /// Set the level for this node
  void set_level(int level);

  /// Return the indices of the level-0 node containing this node
  void array (int * iax, int *iay, int *iaz) const;

  /// Accumulate array part of an index
  void set_array(int iax, int iay, int iaz);

  /// Return the packed tree bits for each axis
  void tree (int * bx = 0, int *by = 0, int *bz = 0,
             int level=INDEX_UNDEFINED_LEVEL) const;

  /// child index of this node in parent
  void child (int level, int * icx, int * icy, int * icz,
              int min_level = 0) const;

  /// Set the child indicies of this node in the parent
  void set_child(int level, int icx, int icy=0, int icz=0,
                 int min_level = 0);

  /// Set this Index to be the given child of the index
  inline void push_child(int icx, int icy=0, int icz=0,
                  int min_level = 0)
  {
    const int level = this->level();
    set_level (level+1);
    set_child (level+1,icx,icy,icz,min_level);
  }

  void print (std::string msg, int level) const;

  void print (const char * msg,
	      int max_level,
	      int rank,
	      const int nb3[3],
	      bool no_nl,
	      void * simulation = 0) const;

  void write (int ip,
	      const char * msg,
	      int max_level,
	      int rank,
	      const int nb3[3]) const;

  std::string bit_string (int max_level,int rank, const int nb3[3]) const;

  /// Comparison operator required for Charm++ pup()
  friend bool operator < (const Index & x, const Index & y) {
    Index a = x;
    Index b = y;
    a.clean_();
    b.clean_();
    if (a.v_[2] < b.v_[2]) return true;
    if (a.v_[2] > b.v_[2]) return false;
    if (a.v_[1] < b.v_[1]) return true;
    if (a.v_[1] > b.v_[1]) return false;
    return  (a.v_[0] < b.v_[0]);
  }

  ///--------------------
  /// PACKING / UNPACKING
  ///--------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  char * load_data (char * buffer);

private: // methods

  /// Clear tree bits that are associated with levels higher than
  /// the actual level
  void clean_ ();

  int num_bits_(int value) const;
	
  void print_ (FILE * fp,
	       const char * msg,
	       int max_level,
	       int rank,
	       const int nb3[3],
	       bool no_nl) const;

private: // attributes

    union {
      NodeBits a_[3];
      int v_[3];
    };
};

#ifndef TEST
  PUPbytes(Index)
#endif

#ifndef TEST
// public:
//   void pup(PUP::er &p) {
//   }
PUPbytes(NodeBits)
#endif

//----------------------------------------------------------------------
#ifndef TEST
class CkArrayIndexIndex:public CkArrayIndex {
  Index * index_;
public:
  CkArrayIndexIndex(const Index &in)
  {
    index_ = new (index) Index(in);
    nInts=sizeof(Index)/sizeof(int);
  }
};
#endif

#endif /* INDEX_HPP */
