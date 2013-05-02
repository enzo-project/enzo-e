#ifndef INDEX_HPP
#define INDEX_HPP

// ARRAY_BITS + TREE_BITS + LEVEL_AXIS_BITS == 32

#define INDEX_MAX_ARRAY_BITS       10
#define INDEX_MAX_ARRAY_INDEX    1024 /* 2 ** INDEX_MAX_ARRAY_BITS */
#define INDEX_MAX_TREE_BITS        20
#define INDEX_MAX_LEVEL_AXIS_BITS   2
#define INDEX_MAX_LEVEL_AXIS_RANGE  4 /* 2 ** AXIS_BITS2 */

struct BIndex {
  unsigned array : INDEX_MAX_ARRAY_BITS; // maximum 1024 x 1024 x 1024 root blocks
  unsigned  tree : INDEX_MAX_TREE_BITS; // maximum INDEX_MAX_TREE_BITS levels / bits
  unsigned level : INDEX_MAX_LEVEL_AXIS_BITS; // maximum 32 level specification
#ifdef CONFIG_USE_CHARM
#ifndef TEST
public:
  void pup(PUP::er &p) {
    
  }
#endif
#endif /* CONFIG_USE_CHARM */
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

  Index() { clear(); }

  Index(int ix, int iy, int iz) 
  {
    clear();
    set_array(ix,iy,iz);
  }

#ifdef CONFIG_USE_CHARM
#ifndef TEST
  void pup(PUP::er &p) {
    p(v_,3);
  }
#endif
#endif /* CONFIG_USE_CHARM */

  void clear () 
  {
    for (int axis=0; axis<3; axis++) {
      a_[axis].array = 0;
      a_[axis].tree = 0;
      a_[axis].level = 0;
    }
  }

  void child (int il, int * icx, int * icy, int * icz) const
  {
    if (icx) (*icx) = (a_[0].tree >> il) & 1;
    if (icy) (*icy) = (a_[1].tree >> il) & 1;
    if (icz) (*icz) = (a_[2].tree >> il) & 1;
  }

  bool is_root() const
  { return (v_[0]==0 && v_[1]==0 && v_[2]==0); }

  void array (int * ix, int *iy, int *iz) const
  { 
    if (ix) (*ix) = a_[0].array;
    if (iy) (*iy) = a_[1].array;
    if (iz) (*iz) = a_[2].array;
  }

  int level () const
  { return a_[0].level + INDEX_MAX_LEVEL_AXIS_RANGE*
      (    a_[1].level + INDEX_MAX_LEVEL_AXIS_RANGE*
	   a_[2].level );  }

  unsigned value (int i) const
  { return v_[i];  }

  void set_level(int il)
  { 
    a_[0].level = (il>>0*INDEX_MAX_LEVEL_AXIS_BITS) & 3;
    a_[1].level = (il>>1*INDEX_MAX_LEVEL_AXIS_BITS) & 3;
    a_[2].level = (il>>2*INDEX_MAX_LEVEL_AXIS_BITS) & 3;
    clean();
  }

  /// Clear tree bits that are associated with levels higher than
  /// the actual level
  void clean ()
  {
    for (int il = level()+1; il<INDEX_MAX_TREE_BITS; il++) {
      set_tree(il,0,0,0);
    }
  }

  /// Accumulate array part of an index
  void set_array(int ix, int iy, int iz)
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

  int tree (int axis) const
  { return a_[axis].tree; }
  
  /// Accumulate a level of the tree
  void set_tree(int il, int ix, int iy=0, int iz=0)
  {
    ASSERT ("Index::set_tree","level must be at least 1",il>0);
    --il;
    a_[0].tree = (a_[0].tree & ~(1<<il)) | (ix<<il);
    a_[1].tree = (a_[1].tree & ~(1<<il)) | (iy<<il);
    a_[2].tree = (a_[2].tree & ~(1<<il)) | (iz<<il);
  }

  void print (const char * msg = "\0") const
  {
    printf ("INDEX %s: ", msg);
    printf ("L [ %d ] ",a_[0].level + INDEX_MAX_LEVEL_AXIS_RANGE*
	    (           a_[1].level + INDEX_MAX_LEVEL_AXIS_RANGE*
			a_[2].level ));
    printf ("T [ ");
    for (int axis=0; axis<3; axis++) {
    printf ("%x ",a_[axis].tree);
    }
    printf ("] ");
    printf ("A [ ");
    for (int axis=0; axis<3; axis++) {
    printf ("%d ",a_[axis].array);
    }
    printf ("] ");
    printf ("[%08X-%08X-%08X]",v_[0],v_[1],v_[2]);
    printf ("\n");
  }
private:
    union {
      BIndex a_[3];
      unsigned v_[3];
    };
};

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
    index[0] = in.value(0);
    index[1] = in.value(1);
    index[2] = in.value(2);
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  }
  //Not required, but convenient: cast-to-foo operators
  operator Index &() {return index_;}
  operator const Index &() const {return index_;}
};
#endif
#endif /* INDEX_HPP */

#endif /* CONFIG_USE_CHARM */
