// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_BlockTrace.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-04-21
/// @brief    [\ref Mesh] Declaration of the BlockTrace class
///
/// This class is used by MethodOutput to help in performing a
/// sequential depth-first traversal of a section of a parallel
/// distributed array of octrees.  It stores the "trace" (sequence of
/// ancestors) of Index's from the root block to the current Block in
/// the traversal, provides a method for updating the trace for the
/// next() Block in the octree, and stores the tree's root block
/// (index_root) and "home" block (the assigned block performing the
/// writing in MethodOutput)

#ifndef MESH_BLOCK_TRACE_HPP
#define MESH_BLOCK_TRACE_HPP

class BlockTrace {

  /// @class    BlockTrace
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh]

public: // interface

  /// Default Constructor
  BlockTrace() throw()
  {}

  /// Constructor for BlockTrace object with given root Index
  BlockTrace(int rank, int index_min[3], int index_max[3]) throw()
    : index_home_(),
      index_root_(),
      index_min_(),
      index_max_(),
      index_curr_(),
      child_stack_(),
      index_stack_()
  {
    for (int i=0; i<3; i++) {
      index_min_[i] = (i<rank) ? index_min[i] : 0;
      index_max_[i] = (i<rank) ? index_max[i] : 1;
      index_curr_[i] = index_min_[i];
    }
    index_home_.set_array(index_min[0],index_min[1],index_min[2]);
    index_home_.set_level(0);
    index_root_ = index_home_;
    ASSERT1("BlockTrace::BlockTrace()",
            "rank must satisfy 1 <= rank = %d <= 3",
            rank,
            (1 <= rank && rank <= 3));
    num_children_ = (rank == 1) ? 2 : (rank == 2) ? 4 : 8;
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    p | num_children_;
    p | index_home_;
    p | index_root_;
    p | child_stack_;
    p | index_stack_;
  };

  /// Go to the "next" Index in the tree, given whether the current
  /// index is a leaf or not.  Return false when we are done traversing
  /// the tree
  bool next (bool is_leaf)
  {
    if (is_leaf) {
      // If this is a leaf node
      while ((child_stack_.size() > 0) &&
             child_stack_.back() >= (num_children_ - 1)) {
        // ... back up to first incompletely visited node or root
        child_stack_.pop_back();
        index_stack_.pop_back();
      }
      // ASSERT child_stack_.size() == 0
      // or child_stack_.back() < num_children - 1
      if (child_stack_.size() > 0) {
        // ... update node child and index
        int ic = ++child_stack_.back();
        int ix = (ic & 1) ? 1 : 0;
        int iy = (ic & 2) ? 1 : 0;
        int iz = (ic & 4) ? 1 : 0;
        // ... get next Child index
        const int level = index_stack_.back().level();
        index_stack_.back().set_child(level,ix,iy,iz);

      } else {
        // next index in the array, if any
        // reset index_curr_ = index_min_ if done
        if (++index_curr_[0] >= index_max_[0]) {
          index_curr_[0] = index_min_[0];
          if (++index_curr_[1] >= index_max_[1]) {
            index_curr_[1] = index_min_[1];
            if (++index_curr_[2] >= index_max_[2]) {
              index_curr_[2] = index_min_[2];
            }
          }
        }
        // update root index to current
        index_root_.set_array
          (index_curr_[0], index_curr_[1], index_curr_[2]);
      }
    } else {
      // push back 0 if not a leaf
      const int level = index_stack_.size();
      Index index = top();
      index.set_level(level+1);
      index.set_child(level+1,0,0,0);

      child_stack_.push_back(0);
      index_stack_.push_back(index);
    }
    return ((child_stack_.size() == 0) &&
            (index_curr_[0]==index_min_[0]) &&
            (index_curr_[1]==index_min_[1]) &&
            (index_curr_[2]==index_min_[2]));
  }

  /// Return the home Index of the BlockTrace
  Index home ()
  { return index_home_; }

  /// Return the top (current) Index of the BlockTrace
  Index top ()
  { return (index_stack_.size() > 0) ? index_stack_.back() : index_root_; }

  bool is_empty()
  { return index_stack_.size() > 0; }

  void print (const char * message)
  {
    CkPrintf ("DEBUG_PRINT BlockTrace %s\n",message);
    int level;
    int a3[3] = {0};
    level = index_home_.level();
    index_home_.array(a3,a3+1,a3+2);
    CkPrintf ("DEBUG_PRINT home level %d array(%d %d %d) \n",
              level,a3[0],a3[1],a3[2]);
    level = index_root_.level();
    index_root_.array(a3,a3+1,a3+2);
    CkPrintf ("DEBUG_PRINT root level %d array(%d %d %d) \n",
              level,a3[0],a3[1],a3[2]);
    for (size_t i=0; i<index_stack_.size(); i++) {
      int c3[3];
      Index index = index_stack_[i];
      int child = child_stack_[i];
      level = index.level();
      a3[0]=a3[1]=a3[2]=0;
      index.array(a3,a3+1,a3+2);
      c3[0]=c3[1]=c3[2]=0;
      index.child(level,c3,c3+1,c3+2);
      CkPrintf ("            %lu level %d a(%d %d %d) c(%d %d %d) [%d]\n",
                i,level,a3[0],a3[1],a3[2],c3[0],c3[1],c3[2],child);
    }
  }

  /// Return the number of bytes required to serialize the data object
  int data_size () const;

  /// Serialize the object into the provided empty memory buffer.
  char * save_data (char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  char * load_data (char * buffer);

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Number of child Blocks for each Block, computed from rank parameter
  int num_children_;

  /// Home Index of the array-of-trees section being traversed
  Index index_home_;

  /// Root Index of the tree being traversed
  Index index_root_;
  
  /// Range and current array indices
  int index_min_[3];
  int index_max_[3];
  int index_curr_[3];

  /// Stack of last visited child indices
  std::vector <int>   child_stack_;

  /// Stack of Block Indexes
  std::vector <Index> index_stack_;

};

#endif /* MESH_BLOCK_TRACE_HPP */

