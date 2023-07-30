// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodOrder.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-01
/// @brief    [\ref Problem] Declaration of the MethodOrder class for
///           generating a space-filling curve (SFC) ordering of blocks in the hierarchy

#ifndef PROBLEM_METHOD_ORDER_HPP
#define PROBLEM_METHOD_ORDER_HPP

class MethodOrder : public Method {

  /// @class    MethodOrder
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  MethodOrder(int min_level) throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(MethodOrder);
  
  /// Charm++ PUP::able migration constructor
  MethodOrder (CkMigrateMessage *m)
    : Method (m)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    Method::pup(p);
    p | type_;
    p | is_index_;
    p | is_count_;
    p | is_next_;
    p | is_windex_;
    p | is_wcount_;
    p | is_count_child_;
    p | is_wcount_child_;
    p | is_sync_index_;
    p | is_sync_count_;
    p | min_level_;
  }

public: // virtual methods

  /// Apply the method to determine the ordering of blocks
  virtual void compute( Block * block) throw();

  virtual std::string name () throw ()
  { return "order"; }

public: // methods

  // accumulate block counts and weights
  void accum_count(Block * block, int count, double wcount, const int ic3[3],
                   int sync_stop = 0);
  // accumulate the final block index in the ordering
  void accum_index(Block * block, int index, int count, double windex, double wcount,
                   int sync_stop = 0);

private: // methods

  /// Create and initialize a MsgOrder object
  MsgOrder * create_msg_order_(Block * block);

  void compute_complete_(Block * block);

  /// The block's ordering index
  long long & index_(Block * block);

  /// The count of participating blocks
  long long & count_(Block * block);

  /// The weighted index of the block
  double & windex_(Block * block);

  /// The weighted count of the block
  double & wcount_(Block * block);

  /// Return the given Block's child block count
  long long & count_child_(Block * block, int index);

  /// Return the given Block's child weight
  double & wcount_child_(Block * block, int index);

  /// Return the next block in the ordering
  Index & next_(Block * block);

  /// Return the Block's ordering index
  Sync & sync_index_(Block * block);

  /// Return the Block's weight (including self)
  Sync & sync_count_(Block * block);

private: // functions

  inline int child_order_(const int ic3[3])
  { return ic3[0] + 2*(ic3[1] + 2*ic3[2]);  }

  inline void child_order_(int ic3[3], int i)
  { int j=i;
    ic3[0] = j & 1;    j = (j >> 1);
    ic3[1] = j & 1;    j = (j >> 1);
    ic3[2] = j & 1;
    ASSERT ("child_order_(ic3)",
            "child_order_(ic3,i) inconsistent with child_order_(ic3)",
            (i == child_order_(ic3)));
  }

  void set_next_ (Block * block);

  inline double weight_(Block * block)
  { return 0.5; }

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// enum type of ordering
  enum class Type { morton };
  Type type_;

  /// Block Scalar<long long> index
  int is_index_;
  /// Block Scalar<long long> count
  int is_count_;
  /// Block Scalar<double> weighted index
  int is_windex_;
  /// Block Scalar<double> weighted count
  int is_wcount_;
  /// Block Scalar<int> child weight (array of size cello::num_children())
  int is_count_child_;
  /// Block Scalar<double> child weight (array of size cello::num_children())
  int is_wcount_child_;
  /// Block Scalar<Index> next Index in ordering
  int is_next_;
  /// Block Scalar<Sync> sync counter for index (coarse->fine)
  int is_sync_index_;
  /// Block Scalar<sync> sync counter for count (fine->coarse)
  int is_sync_count_;

  /// Minimum refinement level for ordering; may be < 0
  int min_level_;
};

#endif /* PROBLEM_METHOD_ORDER_HPP */

